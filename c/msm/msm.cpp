#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <memory.h>
#include <iostream>
#include <vector>
#include <numeric>
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

    // finish processing remaining elements
    finalize();

    // combine all the results into a single point
    reduce(r);
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

    // Separate bucket_ids and point_idxs from batch_buckets_and_points
    std::vector<uint32_t> bucket_ids, point_idxs;
    for (const auto &bp : batchBucketsAndPoints) {
        bucket_ids.push_back(bp.first);
        point_idxs.push_back(bp.second);
    }

    // Perform batch addition
    batchAdder.batch_add_indexed(buckets, bucket_ids, curPoints, point_idxs);

    // Clean up current batch
    bitmap.clear();
    batchBucketsAndPoints.clear();

    // Memorize the last point which is the current processing point
    auto slicing_point = curPoints.back();
    curPoints.pop_back();  // Remove the last point
    curPoints.clear();

    // Process collisions
    int next_pos = 0;
    for (int i = 0; i < collisionBucketsAndPoints.size(); i++) {
        auto bucket_and_point = collisionBucketsAndPoints[i];
        auto bucket_id = bucket_and_point.first; // or bucket_and_point.get<0>(), etc.
        auto point = bucket_and_point.second; // or bucket_and_point.get<1>(), etc.


        if (bitmap.testAndSet(bucket_id)) {
            // collision found
            std::swap(collisionBucketsAndPoints[next_pos], collisionBucketsAndPoints[i]);
            next_pos += 1;
        } else {
            collisionBucketsAndPoints.push_back(std::make_pair(bucket_id, curPoints.size()));
            curPoints.push_back(point);
        }
    }


    // Remove processed collisions
    collisionBucketsAndPoints.resize(next_pos);

    // Push back the slicing_point to curPoints
    curPoints.push_back(slicing_point);
}

template <typename Curve>
void MSM<Curve>::finalize() {
    processBatch();
    while ( !( (batchBucketsAndPoints.size() == 0)  && (collisionBucketsAndPoints.size() == 0) ) )  {
        processBatch();
    }
}

template <typename Curve>
void MSM<Curve>::reduce(typename Curve::Point &r) {
    std::vector<uint32_t> window_starts(nSlices);
    std::iota(window_starts.begin(), window_starts.end(), 0);

    std::vector<typename Curve::Point> window_sums;
    window_sums.reserve(window_starts.size());
    for (auto &w_start : window_starts) {
        size_t bucket_start = static_cast<size_t>(w_start << bitsPerSlice);
        size_t bucket_end = bucket_start + (1 << bitsPerSlice);
        window_sums.push_back(innerWindowReduce(bucket_start, bucket_end));
    }


    r = intraWindowReduce(window_sums);
}

template <typename Curve>
typename Curve::Point MSM<Curve>::intraWindowReduce(std::vector<typename Curve::Point> window_sums) {
    typename Curve::Point lowest = window_sums.front();
    typename Curve::Point total = Curve::zero();

    for (auto it = window_sums.rbegin(); it != window_sums.rend() - 1; ++it) {
        total += *it;
        for (int i = 0; i < bitsPerSlice; ++i) {
            total.double_in_place();
        }
    }

    return lowest + total;
}



template <typename Curve>
typename Curve::Point MSM<Curve>::innerWindowReduce(size_t start, size_t end) {
    typename Curve::Point running_sum = Curve::zero();
    typename Curve::Point res = Curve::zero();

    for (auto it = buckets.rbegin() + (buckets.size() - end); it != buckets.rbegin() + (buckets.size() - start); ++it) {
        running_sum.add_assign_mixed(*it);
        res += running_sum;
    }
    return res;
}

