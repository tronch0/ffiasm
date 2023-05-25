#ifndef PAR_MULTIEXP2
#define PAR_MULTIEXP2

#include "bitmap.h"
#include "batch_adder.h"

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
    Curve &g;
    uint32_t nSlices;

    std::vector<typename Curve::Point> curPoints;

    std::vector<std::pair<uint32_t, uint32_t>> batchBucketsAndPoints;
    uint32_t maxBatchCnt;
    std::vector<std::pair<uint32_t, typename Curve::Point>> collisionBucketsAndPoints;
    uint32_t maxCollisionCnt;

    Bitmap bitmap;
    std::vector<typename Curve::Point> buckets;

    BatchAdder<Curve, typename Curve::PointAffine> batchAdder;

public:

    MSM(Curve &_g): g(_g) { batchAdder = BatchAdder<Curve, typename Curve::PointAffine>();}

    void run(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads= 0);
    std::vector<uint32_t> slice(uint32_t scalarIdx, uint32_t scalarSize, uint32_t bitsPerChunk);
    void processPointAndSlices(uint32_t baseIdx, std::vector<uint32_t>& slices);
    void processSlices(uint32_t bucket_id, uint32_t currentPointIdx) ;
    void processBatch();
    void finalize();
    void reduce(typename Curve::Point &r);
    typename Curve::Point innerWindowReduce(size_t start, size_t end);
    typename Curve::Point intraWindowReduce(std::vector<typename Curve::Point> window_sums);

    };

#include "msm.cpp"

#endif // PAR_MULTIEXP2
