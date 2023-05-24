#ifndef PAR_MULTIEXP2
#define PAR_MULTIEXP2

#include "bitmap.h"

#define PME2_PACK_FACTOR 2
#define PME2_MAX_CHUNK_SIZE_BITS 16
#define PME2_MIN_CHUNK_SIZE_BITS 2

template <typename Curve>
class MSM {

    typename Curve::PointAffine *bases;
    uint8_t* scalars;
    uint32_t scalarSize;
    uint32_t n;
    uint32_t nThreads;
    uint32_t bitsPerSlice;
    uint32_t nSlices;
    Curve &g;

    uint32_t cur_batch_cnt;
    uint32_t max_batch_cnt;

    Bitmap bitmap;

    uint32_t collision_cnt;
    uint32_t max_collision_cnt;

    std::vector<typename Curve::Point> cur_points;
    uint32_t cur_points_cnt;



    std::vector<typename Curve::Point> buckets;

    BatchAdder<Curve, typename Curve::PointAffine> ba;

    // New Body
    // ***************

//    num_windows: u32,
//            window_bits: u32,
//            bucket_bits: u32,
//            max_batch_cnt: u32,          // max slices allowed in a batch
//    max_collision_cnt: u32,
//            buckets: Vec<GroupAffine<P>>, // size (num_windows << window_bits) * 2
//
//    // current batch state
//    bitmap: Bitmap,
//            batch_buckets_and_points: Vec<(u32, u32)>,
//            collision_buckets_and_points: Vec<(u32, GroupAffine<P>)>,
//            cur_points: Vec<GroupAffine<P>>, // points of current batch, size batch_size
//
//    // batch affine adder
//    batch_adder: BatchAdder<P>,

public:
    // ctor:

//    pub fn new(
//    scalar_bits: u32,
//            window_bits: u32,
//            max_batch_cnt: u32,     // default: 4096
//    max_collision_cnt: u32, // default: 128
//    ) -> BucketMSM<P> {
//        let num_windows = (scalar_bits + window_bits - 1) / window_bits;
//        let batch_size = std::cmp::max(8192, max_batch_cnt);
//        let bucket_bits = window_bits - 1; // half buckets needed because of signed-bucket-index
//        let bucket_size = num_windows << bucket_bits;
//        // size of batch_adder will be the max of batch_size and num_windows * groups per window
//        let batch_adder_size = std::cmp::max(batch_size, bucket_size >> GROUP_SIZE_IN_BITS);
//
//        BucketMSM {
//                num_windows,
//                window_bits,
//                bucket_bits,
//                max_batch_cnt,
//                max_collision_cnt,
//                buckets: vec![GroupAffine::<P>::zero(); bucket_size as usize],
//
//                bitmap: Bitmap::new(bucket_size as usize / 32),
//                batch_buckets_and_points: Vec::with_capacity(batch_size as usize),
//                collision_buckets_and_points: Vec::with_capacity(max_collision_cnt as usize),
//                cur_points: vec![GroupAffine::<P>::zero(); batch_size  as usize],
//
//                batch_adder: BatchAdder::new(batch_adder_size as usize),
//        }
//    }





    MSM(Curve &_g): g(_g) {ba = BatchAdder<Curve, typename Curve::PointAffine>();}

    void run(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads= 0);
    std::vector<uint32_t> slice(uint32_t scalarIdx, uint32_t scalarSize, uint32_t bitsPerChunk);
    void process_batch();
};

#include "msm.cpp"

#endif // PAR_MULTIEXP2
