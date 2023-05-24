#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <memory.h>
#include <iostream>
#include <vector>
#include "../misc.hpp"
#include "msm.hpp"

template <typename Curve>
void MSM<Curve>::run(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads) {

    // calculate and stored all members of the MSM.
    // *******************************************

    // double check if to remove any constraints from here (all consts are from previous implementation)
    bitsPerSlice = log2(n / PME2_PACK_FACTOR);
    if (bitsPerSlice > PME2_MAX_CHUNK_SIZE_BITS) bitsPerSlice = PME2_MAX_CHUNK_SIZE_BITS;
    if (bitsPerSlice < PME2_MIN_CHUNK_SIZE_BITS) bitsPerSlice = PME2_MIN_CHUNK_SIZE_BITS;
    nSlices = ((scalarSize * 8 - 1 ) / bitsPerSlice) + 1;

    bases = _bases;
    scalars = _scalars;
    scalarSize = _scalarSize;
    n = _n;
    nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;


    for(uint32_t i=0; i<n; i++) {           // the for is iterating over the base points
        std::vector<uint32_t> slices = slice(i, scalarSize, bitsPerSlice);
        processPointAndSlices(i, slices);
    }

    // process complete

    // reduce
}

template <typename Curve>
std::vector<uint32_t> MSM<Curve>::slice(uint32_t scalarIdx, uint32_t scalarSize, uint32_t bitsPerChunk) {
    std::vector<uint32_t> slices;
    uint32_t chunksCount = (scalarSize * 8 + bitsPerChunk - 1) / bitsPerChunk; // Calculate chunks count, same as slices in rust
    slices.resize(chunksCount);

    for (uint32_t chunkIdx = 0; chunkIdx < chunksCount; ++chunkIdx) {
        uint32_t bitStart = chunkIdx * bitsPerChunk;
        uint32_t byteStart = bitStart / 8;
        uint32_t efectiveBitsPerChunk = bitsPerChunk;
        if (byteStart > scalarSize-8) byteStart = scalarSize - 8;
        if (bitStart + bitsPerChunk > scalarSize*8) efectiveBitsPerChunk = scalarSize*8 - bitStart;
        uint32_t shift = bitStart - byteStart*8;
        uint64_t v = *(uint64_t *)(scalars + scalarIdx*scalarSize + byteStart);
        v = v >> shift;
        v = v & ( (1 << efectiveBitsPerChunk) - 1);
        slices[chunkIdx] = uint32_t(v);
    }
    return slices;
}

template <typename Curve>
void MSM<Curve>::processPointAndSlices(uint32_t baseIdx, std::vector<uint32_t>& slices) {
    assert(nSlices == slices.size() && "slices size should equal num_windows");  // TODO Remove later for efficiency

    typename Curve::PointAffine point = bases[baseIdx];

    curPoints.push_back(point);

    for (int win = 0; win < slices.size(); win++) {
        if (slices[win] > 0) {
            uint32_t bucket_id = (win << bitsPerSlice) + slices[win] - 1; // skip slice == 0
            processSlices(bucket_id, curPoints.size() - 1);
        }
    }
}

template <typename Curve>
void MSM<Curve>::processSlices(uint32_t bucket_id, uint32_t currentPointIdx) {
    // Check if the bucket_id bit is not set
    if (!bitmap.testAndSet(bucket_id)) {
        // If no collision found, add point to current batch
        batchBucketsAndPoints.push_back(std::make_pair(bucket_id, currentPointIdx));
    } else {
        // In case of collision, add it to the collision list
        collisionBucketsAndPoints.push_back(std::make_pair(bucket_id, curPoints[currentPointIdx]));
    }

    // If the count of collisions or the batch size reach their maximum limits, process the batch
    if (collisionBucketsAndPoints.size() >= maxCollisionCnt ||
            batchBucketsAndPoints.size() >= maxBatchCnt) {
        processBatch();
    }
}

template <typename Curve>
void MSM<Curve>::processBatch(){
    if (batchBucketsAndPoints.empty()) {
        return;
    }

    
}






//template <typename Curve>
//void MSM<Curve>::process_batch() {
//
//}

//template <typename Curve>
//void MSM<Curve>::reduce() {
//
//}


