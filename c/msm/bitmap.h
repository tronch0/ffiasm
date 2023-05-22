#ifndef RAPIDSNARK_BITMAP_H
#define RAPIDSNARK_BITMAP_H

#include <vector>

class Bitmap {
private:
    size_t size;
    std::vector<uint64_t> data;

public:
    Bitmap(size_t size) : size(size), data((size + 63) / 64, 0) {}

    bool testAndSet(uint32_t bucket) {
        if (bucket >= size) {
            // Handle out-of-range bucket
            return false;
        }

        uint32_t block = bucket / 64;
        uint32_t offset = bucket % 64;
        uint64_t bitMask = 1ULL << offset;

        bool current = (data[block] & bitMask) != 0;
        data[block] |= bitMask;
        return current;
    }

    void clear() {
        data.assign((size + 63) / 64, 0);
    }
};
#endif //RAPIDSNARK_BITMAP_H


//#include <vector>
//#include <bitset>

//class Bitmap {
//private:
//    size_t size;
//    std::vector<std::bitset<32>> data;
//
//public:
//    Bitmap(size_t size) : size(size), data((size + 31) / 32) {}
//
//    bool testAndSet(uint32_t bucket) {
//        size_t word = bucket / 32;
//        size_t bit = bucket % 32;
//
//        if (data[word].test(bit)) {
//            return true;
//        }
//
//        data[word].set(bit);
//        return false;
//    }
//
//    void clear() {
//        for (auto& word : data) {
//            word.reset();
//        }
//    }
//};


