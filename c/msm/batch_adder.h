#ifndef RAPIDSNARK_BATCH_ADDER_H
#define RAPIDSNARK_BATCH_ADDER_H
#include <vector>
#include <cassert>
#include <array>

template <typename Curve, typename Field, size_t MaxBatchSize>
class BatchAdder {

public:
    Curve &g;

    std::array<typename Field::Element, MaxBatchSize> inverses; //    Alternative: std::vector<typename Field::Element> inverses;
//    size_t batch_size;
    typename Field::Element inverse_state;

    BatchAdder(Curve& aCurve)
            : g(aCurve), inverse_state(aCurve.F.one()) {
        size_t max_batch_cnt = std::extent<decltype(inverses)>::value;
//        batch_size = max_batch_cnt;
        //  #pragma omp parallel for    TODO:Uncomment this line to parallelize the for loop
        for (int i=0; i<max_batch_cnt; i++) {
            g.F.copy(inverses[i], g.F.one());
        }
    }

    void batch_add_indexed(
//            typename Curve::PointAffine* dest,
//            const typename Curve::PointAffine* src,
//
//            const std::vector<uint32_t>  &dest_indices,
//            const std::vector<uint32_t>  &src_indices

            std::vector<typename Curve::PointAffine>& dest,
            std::vector<typename Curve::PointAffine>& src,
            std::vector<uint32_t>& dest_indices,
            std::vector<uint32_t>& src_indices

    ) {
        size_t destLen = dest.size();
        size_t srcLen = src.size();
        size_t destIdxLen = dest_indices.size();
        size_t srcIdxLen = src_indices.size();

//        size_t inversesLen = sizeof(inverses) / sizeof(inverses[0]);
//        size_t destLen = sizeof(dest) / sizeof(dest[0]);
//        size_t srcLen = sizeof(src) / sizeof(src[0]);
//
//        size_t destIdxLen = dest_indices.size();
//        size_t srcIdxLen = src_indices.size();

        assert(destLen == srcLen && "insufficient entries in dest array");
        assert(destLen <= inverses.size() && "input length exceeds the max_batch_cnt, please increase max_batch_cnt during initialization!");
        assert(destIdxLen == srcIdxLen && "length of dest_indexes and src_indexes don't match!");

        reset();
        for (size_t i = 0; i < destIdxLen; ++i) {
            size_t dest_idx = dest_indices[i];
            size_t src_idx = src_indices[i];
            batch_add_phase_one(dest[dest_idx], src[src_idx], i);
        }
        inverse();
        for (size_t i = destIdxLen - 1; i >= 0; --i) {
            size_t dest_idx = dest_indices[i];
            size_t src_idx = src_indices[i];
            batch_add_phase_two(dest[dest_idx], src[src_idx], i);
        }
    }

    void inverse() {
//        inverse_state = inverse_state.inverse().unwrap();
        g.F.inv(inverse_state, inverse_state);

    }

    void reset() {
        g.F.copy(inverse_state, g.F.one());
    }

private:

    /// Two-pass batch affine addition
    ///   - 1st pass calculates from left to right
    ///      - inverse_state: accumulated product of deltaX
    ///      - inverses[]: accumulated product left to a point
    ///   - call inverse()
    ///   - 2nd pass calculates from right to left
    ///      - slope s and ss from state
    ///      - inverse_state = inverse_state * deltaX
    ///      - addition result acc

    void batch_add_phase_one(typename Curve::PointAffine &p, typename Curve::PointAffine &q, size_t idx) {
        assert(idx < inverses.size() && "index exceeds the max_batch_cnt, please increase max_batch_cnt during initialization!");

        if (g.isZero(p) || g.isZero(q)) {
            return;
        }

        typename Field::Element delta_x;
        typename Field::Element delta_y;

        g.F.sub(delta_x, q.x, p.x);
        if (g.F.isZero(delta_x)) {
            g.F.sub(delta_y, q.y, p.y);
            if (!g.F.isZero(delta_y)) { // p = -q, return
                return;
            }

            // if P == Q
            // if delta_x is zero, we need to invert 2y
            g.F.add(delta_x, q.y, p.y); // q.Y + q.Y;
        }

        if (g.F.isZero(inverse_state)) {
            g.F.copy(inverses[idx], g.F.one());
            g.F.copy(inverse_state, delta_x);
        } else {
            g.F.copy(inverses[idx], inverse_state);
            g.F.mul(inverse_state, inverse_state, delta_x);
        }
    }

    /// should call inverse() between phase_one and phase_two
    void batch_add_phase_two(typename Curve::PointAffine& p, typename Curve::PointAffine& q, size_t idx) {
        assert(idx < inverses.size() && "Size mismatch between dest and indices!");

        if (g.isZero(p) || g.isZero(q)) {
            if (!g.isZero(q)) {
                g.copy(p, q);
            }
            return;
        }

        typename Field::Element inverse;
        g.F.copy(inverse, inverses[idx]);

        g.F.mul(inverse, inverse, inverse_state);

        typename Field::Element delta_x;
        typename Field::Element delta_y;

        g.F.sub(delta_x, q.x, p.x);  //q.X - p.X;
        g.F.sub(delta_y, q.y, p.y);  //q.Y - p.Y;

        if (g.F.isZero(delta_x)) {
            if (!g.F.isZero(delta_y)) {
                g.copy(p, g.zeroAffine());
                return;
            }

            g.F.square(delta_y, q.x);
            g.F.add(delta_y, delta_y, delta_y); // delta_y += delta_y;
            g.F.add(delta_y, delta_y, delta_y); // delta_y += delta_y;

            g.F.square(delta_x, q.y); // delta_x = q.Y.doubled();
        }

        g.F.mul(inverse_state, inverse_state, delta_x);

        typename Field::Element s;
        typename Field::Element ss;

        g.F.mul(s, delta_y, inverse);
        g.F.square(ss, s);

        g.F.sub(ss, ss, q.x);
        g.F.sub(ss, ss, p.x);

        g.F.sub(delta_x, q.x, p.x);  //q.X - p.X;

        g.F.sub(p.y, s, delta_x);

        g.F.sub(p.y, p.y, q.y);
    }
};

#endif //RAPIDSNARK_BATCH_ADDER_H
