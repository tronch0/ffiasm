#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <memory.h>
#include "misc.hpp"
#include "msm_improve.hpp"

/*
template<typename G, typename G1Affine, typename BaseField>
class BucketMSM {
    private:
        uint32_t num_windows;
        uint32_t window_bits;

        std::vector<G1Affine> buckets;
        G results;

        std::vector<G1Affine> cur_points;
        uint32_t cur_points_cnt;
        std::vector<uint64_t> cur_bkt_pnt_pairs;
        uint32_t cur_batch_cnt;
        uint32_t max_batch_cnt;

        BaseField inverse_state;
        std::vector<BaseField> inverses;

        Bitmap bitmap;  // bitmap for the current window
        uint32_t collision_cnt;
        uint32_t max_collision_cnt;

        std::vector<uint64_t> collision_list_nodes;
        uint32_t available_list;
        uint32_t processing_list;
        uint32_t unprocessed_list;

    public:
        // constructor and functions
};
*/

// *** Test unit testing framework ***
// ***********************************
int addTwoNumbers() {
    int num1 = 3;
    int num2 = 7;
    return num1 + num2;
}

bool isEvenNumber() {
    int num = 4;
    return (num % 2) == 0;
}


//void msm_slice(const vector<uint8_t>& scalar, vector<uint32_t>& slices, const uint32_t window_bits) {
//    cpp_int temp = 0;
//    for(int i = scalar.size() - 1; i >= 0; --i)
//        temp = temp * 256 + scalar[i];
//
//    for(size_t i = 0; i < slices.size(); ++i) {
//        slices[i] = static_cast<uint32_t>(temp % (1 << window_bits));
//        temp /= (cpp_int(1) << window_bits);
//    }
//}



template <typename Curve>
void MSM<Curve>::initAccs() {
    #pragma omp parallel for
    for (int i=0; i<nThreads*accsPerChunk; i++) {
        g.copy(accs[i].p, g.zero());
    }
}

template <typename Curve>
uint32_t MSM<Curve>::getChunk(uint32_t scalarIdx, uint32_t chunkIdx) {
    uint32_t bitStart = chunkIdx*bitsPerChunk;
    uint32_t byteStart = bitStart/8;
    uint32_t efectiveBitsPerChunk = bitsPerChunk;
    if (byteStart > scalarSize-8) byteStart = scalarSize - 8;
    if (bitStart + bitsPerChunk > scalarSize*8) efectiveBitsPerChunk = scalarSize*8 - bitStart;
    uint32_t shift = bitStart - byteStart*8;
    uint64_t v = *(uint64_t *)(scalars + scalarIdx*scalarSize + byteStart);
    v = v >> shift;
    v = v & ( (1 << efectiveBitsPerChunk) - 1);
    return uint32_t(v);
}

// place into buckets
// go over all the numbers (windowed numbered) in the window/chunk and add them to their corresponding index
template <typename Curve>
void MSM<Curve>::processChunk(uint32_t idChunk) {
    #pragma omp parallel for
    for(uint32_t i=0; i<n; i++) {           // iterate over the window's chunks and add them to their corresponding index
        if (g.isZero(bases[i])) continue;
#ifdef _OPENMP
        int idThread = omp_get_thread_num();
#else
        int idThread = 0;
#endif
        uint32_t chunkValue = getChunk(i, idChunk);
        if (chunkValue) {
            g.add(accs[idThread*accsPerChunk+chunkValue].p, accs[idThread*accsPerChunk+chunkValue].p, bases[i]);
        }
    }
}

// This function takes all chunks and accumulate them to the first chunk's indexes
template <typename Curve>
void MSM<Curve>::packThreads() {
    #pragma omp parallel for
    for(uint32_t i=0; i<accsPerChunk; i++) {
        for(uint32_t j=1; j<nThreads; j++) {
            if (!g.isZero(accs[j*accsPerChunk + i].p)) {
                g.add(accs[i].p, accs[i].p, accs[j*accsPerChunk + i].p);
                g.copy(accs[j*accsPerChunk + i].p, g.zero());
            }
        }
    }
}

template <typename Curve>
void MSM<Curve>::reduce(typename Curve::Point &res, uint32_t nBits) {
    if (nBits==1) {
        g.copy(res, accs[1].p);
        g.copy(accs[1].p, g.zero());
        return;
    }
    uint32_t ndiv2 = 1 << (nBits-1);


    PaddedPoint *sall = new PaddedPoint[nThreads];
    memset(sall, 0, sizeof(PaddedPoint)*nThreads);

    typename Curve::Point p;
    #pragma omp parallel for
    for (uint32_t i = 1; i<ndiv2; i++) {
#ifdef _OPENMP
        int idThread = omp_get_thread_num();
#else
        int idThread = 0;
#endif
        if (!g.isZero(accs[ndiv2 + i].p)) {
            g.add(accs[i].p, accs[i].p, accs[ndiv2 + i].p);
            g.add(sall[idThread].p, sall[idThread].p, accs[ndiv2 + i].p);
            g.copy(accs[ndiv2 + i].p, g.zero());
        }
    }
    for (u_int32_t i=0; i<nThreads; i++) {
        g.add(accs[ndiv2].p, accs[ndiv2].p, sall[i].p);
    }

    typename Curve::Point p1;
    reduce(p1, nBits-1);

    for (u_int32_t i=0; i<nBits-1; i++) g.dbl(accs[ndiv2].p, accs[ndiv2].p);
    g.add(res, p1, accs[ndiv2].p);
    g.copy(accs[ndiv2].p, g.zero());
    delete[] sall;
}

template <typename Curve>
void MSM<Curve>::run(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads) {
#ifdef _OPENMP
    nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
    ThreadLimit threadLimit (nThreads);
#else
    nThreads = 1;
#endif
    bases = _bases;
    scalars = _scalars;
    scalarSize = _scalarSize;
    n = _n;

    if (n==0) {
        g.copy(r, g.zero());
        return;
    }
    if (n==1) {
        g.mulByScalar(r, bases[0], scalars, scalarSize);
        return;
    }
    bitsPerChunk = log2(n / PME2_PACK_FACTOR);
    if (bitsPerChunk > PME2_MAX_CHUNK_SIZE_BITS) bitsPerChunk = PME2_MAX_CHUNK_SIZE_BITS;
    if (bitsPerChunk < PME2_MIN_CHUNK_SIZE_BITS) bitsPerChunk = PME2_MIN_CHUNK_SIZE_BITS;
    nChunks = ((scalarSize*8 - 1 ) / bitsPerChunk)+1;
    accsPerChunk = 1 << bitsPerChunk;  // In the chunks last bit is always zero.

    // prints all information about the process:
    std::cout << "*** MSM Initiate *** (threads: " << nThreads << ", window-bits: " << bitsPerChunk << ", windows: " << nChunks << ", points: " << n << ") \n";

    typename Curve::Point *chunkResults = new typename Curve::Point[nChunks];
    accs = new PaddedPoint[nThreads*accsPerChunk];
    // std::cout << "InitTrees " << "\n";
    initAccs();

    // iterate over the windows/chunks, process all relevant windows for each base point
    for (uint32_t i=0; i<nChunks; i++) {
        // std::cout << "process chunks " << i << "\n";
        processChunk(i);
        // std::cout << "pack " << i << "\n";
        packThreads();
        // std::cout << "reduce " << i << "\n";
        reduce(chunkResults[i], bitsPerChunk);
    }

    delete[] accs;

    g.copy(r, chunkResults[nChunks-1]);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) g.dbl(r,r);
        g.add(r, r, chunkResults[j]);
    }

    delete[] chunkResults;
}
