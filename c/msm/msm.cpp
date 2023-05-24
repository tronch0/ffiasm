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


    for(uint32_t i=0; i<n; i++) {           // the for is iterating over the base points
        std::vector<uint32_t> slices = slice(i, scalarSize, bitsPerSlice);

        // process(i,slices
    }

    // process complete

    // reduce!
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
    typename Curve::PointAffine point = bases[baseIdx];
    cur_points[cur_points_cnt] = point;
    cur_points_cnt++;
    for (int win = 0; win < slices.size(); win++) {
        if (slices[win] > 0) {
            uint32_t bucket_id = (win << bitsPerSlice) + slices[win] - 1;

            if (!bitmap.testAndSet(bucket_id)) {
                const auto& acc = buckets[bucket_id];
                batch_adder::batchAddPhaseOne(acc, point, cur_batch_cnt);

                cur_bkt_pnt_pairs[cur_batch_cnt] = ((uint64_t)bucket_id << 14) + cur_points_cnt - 1;
                cur_batch_cnt++;
            } else {
                int free_node_index = list::pop(collision_list_nodes, available_list);

                assert(free_node_index < 2 * max_collision_cnt);
                collision_list_nodes[free_node_index] = ((uint64_t)bucket_id << 14) + 0x3FFF;
                list::enqueue(collision_list_nodes, unprocessed_list, free_node_index);

                collision_cnt++;
                cur_points[max_batch_cnt + free_node_index] = point;
            }

            if (collision_cnt >= max_collision_cnt || cur_batch_cnt >= max_batch_cnt) {
                process_batch();
            }
        }
    }
}

template <typename Curve>
void MSM<Curve>::process_batch() {

}

template <typename Curve>
void MSM<Curve>::reduce() {

}


