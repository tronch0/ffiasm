#ifndef PAR_MULTIEXP2
#define PAR_MULTIEXP2

#include "bitmap.h"
#include "batch_adder.h"

#define PME2_PACK_FACTOR 2
#define PME2_MAX_CHUNK_SIZE_BITS 16
#define PME2_MIN_CHUNK_SIZE_BITS 2

const size_t MAX_BATCH_SIZE = 1024;
const size_t MAX_COLLISION_SIZE = 2048;

template <typename Curve, typename BaseField>
class MSM {

    typename Curve::PointAffine *bases;
    uint8_t* scalars;
    uint32_t scalarSize;
    uint32_t n;
    uint32_t nThreads;
    uint32_t bitsPerSlice;
    Curve &g;
    uint32_t nSlices;

    std::vector<typename Curve::PointAffine> curPoints;

    std::vector<std::pair<uint32_t, uint32_t>> batchBucketsAndPoints;
    uint32_t maxBatchCnt;
    std::vector<std::pair<uint32_t, typename Curve::PointAffine>> collisionBucketsAndPoints;
    uint32_t maxCollisionCnt;

    Bitmap bitmap;
    std::vector<typename Curve::PointAffine> buckets;

    BatchAdder<Curve, BaseField> batchAdder;

public:

    MSM(Curve &_g): g(_g), bitmap(MAX_BATCH_SIZE), batchAdder(g, MAX_BATCH_SIZE), maxBatchCnt(MAX_BATCH_SIZE), maxCollisionCnt(MAX_COLLISION_SIZE) { }

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
