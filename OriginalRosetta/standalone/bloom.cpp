#include "dst.h"
#include <set>
#include <vector>
#include <cmath>

int main() {
    size_t n=1000000;
    BloomFilter<> bf((size_t)ceil((double)n*8));
    vector<uint32_t> v;
    for (size_t i=0; i<n; ++i) {
        v.push_back(rand());
    }
    clock_t begin, end;
    begin = clock();
    bf.AddKeys(vector<Bitwise>(v.begin(), v.end()));
    end = clock();
    printf("Bloom filter construction took %lfus per key\n", (double)(end-begin)/CLOCKS_PER_SEC*1000000/n);
    set<uint32_t> s(v.begin(), v.end());

    size_t fp = 0, empty = 0;
    for (size_t i=0; i<1000000; ++i) {
        uint32_t key = rand();
        bool full = false, positive = false;

        if (s.find(key) != s.end()) {
            full = true;
        }
        if (bf.Query(key)) {
            positive = true;
        }
        if (not full) {
            ++empty;
            if (positive) {
                ++fp;
            }
        }
        if (full and not positive) {
            printf("FALSE NEGATIVE!!\n");
            return 0;
        }
    }
    printf("FPR: %lf\n", (double)fp/empty);
//    begin = clock();
//    for (size_t i=0; i<1000000; ++i) {
//        uint32_t key = rand();
//
//        if (bf.Query(key)) {
//            ++fp;
//        }
//    }
//    end = clock();
//    printf("printing to avoid optimization: %lu\n", fp);
//    printf("Bloom filter takes %lfus per query\n", (double)(end-begin)/CLOCKS_PER_SEC*1000000/1000000);
}
