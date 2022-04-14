#include <set>
#include <fstream>
#include <iostream>
#include "dst.h"
#ifdef USE_DISK
#include "disk.h"
#endif

#ifndef KEYTYPE
#define KEYTYPE uint64_t
#endif

int main(int argc, char **argv) {
    if (argc < 8) {
        printf("Usage:\n\t%s [keys] [queries-left] [queries-right] [bits-per-key] [SuRF or DST] ([surf hash length] [surf real length]) || ([dfs-diffidence] [bfs-diffidence] [cutoff])\n", argv[0]);
        return 1;
    }

    ifstream key_file, lft_file, rgt_file;
    key_file.open(argv[1]);
    lft_file.open(argv[2]);
    rgt_file.open(argv[3]);
    double bits_per_key = strtod(argv[4], nullptr), dfs_diff=0, bfs_diff=0, cutoff=0, surf_hlen=0, surf_rlen=0;
    bool use_DST = (strcmp(argv[5], "DST") == 0);

    if (use_DST) {
        dfs_diff = strtoull(argv[6], nullptr, 10);
        bfs_diff = strtoull(argv[7], nullptr, 10);
        cutoff = strtoull(argv[8], nullptr, 10);
    }
    else {
        surf_hlen = strtoull(argv[6], nullptr, 10);
        surf_rlen = strtoull(argv[7], nullptr, 10);
    }

#ifndef USE_DISK
    set<KEYTYPE> keyset;
#else
    vector<DiskBlock<KEYTYPE>> files;
    vector<KEYTYPE> keybuffer;
#endif
    vector<Bitwise> keysbwise;
    vector<pair<KEYTYPE, KEYTYPE>> queries;

    size_t maxlen = 0;
    KEYTYPE key, prev_key;
    while (key_file >> key) {
        if (maxlen > 0) {
            assert(prev_key <= key);
        }
        prev_key = key;
#ifndef USE_DISK
        keyset.insert(key);
#else
        keybuffer.push_back(key);
        if (keybuffer.size() == MAX_BUFFER_SIZE) {
            files.emplace_back(keybuffer, "./disk_experiment_binary."+to_string(files.size()), 128);
            keybuffer.clear();
        }
#endif
        keysbwise.push_back(key);
        maxlen = max(maxlen, keysbwise[keysbwise.size()-1].size());
    }
#ifdef USE_DISK
    files.emplace_back(keybuffer, "./disk_experiment_binary."+to_string(files.size()), 128);
    keybuffer.clear();
#endif
    key_file.close();
    size_t nkeys = keysbwise.size();
    KEYTYPE upper;
    vector<pair<KEYTYPE, KEYTYPE>> range_queries;
    vector<size_t> qdist(64, 0);
    for (size_t i=0; i<64; ++i) {
        qdist[i] = 64-i;
    }

    while (lft_file >> key) {
        rgt_file >> upper;
        range_queries.push_back({key, upper});
    }


    if (use_DST) {
        DstFilter<BloomFilter<true>, true> dst_stat(dfs_diff, bfs_diff, [](vector<size_t> x) -> vector<size_t> {for (size_t i=0; i<x.size(); ++i) {x[i]*=1.44;} return x;});
        vector<Bitwise> tmp;
        tmp.push_back(Bitwise(false, maxlen));
        dst_stat.AddKeys(tmp);
        //printf("Done adding fake keys\n");
        for (auto &q: range_queries) {
            if ((size_t)rand()%10000 < 10000 * 10000 / range_queries.size()) {
                (void)dst_stat.Query(q.first, q.second);
            }
        }
        //printf("Making vector...\n");
        qdist = vector<size_t>(dst_stat.qdist_.begin(), dst_stat.qdist_.end());
        printf("Done collecting stats\n");
    }
//    for (size_t i=0; i<64; ++i) {
//        printf("qdist %lu: %lu\n", i, qdist[i]);
//    }

    printf("### BEGIN STATISTICS ###\n");

    DstFilter<BloomFilter<>> *dst=0;
    if (use_DST) {
        dst = new DstFilter<BloomFilter<>>(dfs_diff, bfs_diff, [bits_per_key, nkeys, qdist, cutoff](vector<size_t> x) -> vector<size_t> {return calc_dst(x, bits_per_key, qdist, cutoff);});
    }
    else {
    }
    clock_t begin, end;
    begin = clock();
    if (use_DST) {
        dst->AddKeys(keysbwise);
    }
    else {
    }
    end = clock();
    printf("\tus/Insert:\t%lf\n", (double)1000000*(end-begin)/CLOCKS_PER_SEC/nkeys);
    keysbwise.clear();

    begin = clock();
    size_t empty = 0, fp = 0;
    for (auto &q: range_queries) {
        key = q.first;
        upper = q.second;

        bool full = false, positive = false;
        bool filterAns;
        if (use_DST) {
            filterAns = dst->Query(key, upper);
        }
        else {
        }

        if (filterAns) {
            positive = true;
        }

#ifndef USE_DISK
        auto it = keyset.lower_bound(key);
        if (it != keyset.end() && (*it) < upper) {
            full = true;
        }
        if (not full) {
            ++empty;
            if (positive) {
                ++fp;
            }
        }
        if (full && !positive) {
            cout << "False negative: " << key << " " << upper << " " << *it << "\n";
            assert(false);
        }
#else
        if (filterAns) {
            positive = true;
            for (auto &file: files) {
                if (file.IsRangeFull(key, upper)) {
                    full = true;
                    break;
                }
            }
        }

        if (positive && !full) {
            ++fp;
        }

        if (!full) { // filterAns will also be false
            ++empty;
        }

#endif
    }
    lft_file.close();
    rgt_file.close();
    end = clock();
    printf("\ttotal us/Query:\t%lf\n", (double)1000000*(end-begin)/CLOCKS_PER_SEC/range_queries.size());

    /*
    printf("fpr: %lf\n", (double)fp/empty);
    pair<uint8_t*, size_t> ser = dst.serialize();
    printf("Using %lf bits per key\n", (double)ser.second*8/keysbwise.size());
    */
    printf("\tFPR:\t%lf\n", (double)fp/empty);
    //printf("fn: %lu, fp: %lu, empty: %lu\nfpr: %lf\n", fn, fp, empty, (double)fp/empty);
    pair<uint8_t*, size_t> ser;
    if (use_DST) {
       ser = dst->serialize();
    }
    else {
    }
    printf("\tBPK:\t%lf\n", (double)ser.second*8/nkeys);
    delete[] ser.first;

    volatile bool res=false;
    begin = clock();
    for (auto &q: range_queries) {
        if (use_DST) {
            res = dst->Query(q.first, q.second);
        }
        else {
        }
    }
    end = clock();
    printf("\tus/Query:\t%lf\n", (double)1000000*(end-begin)/CLOCKS_PER_SEC/range_queries.size());

    printf("### END STATISTICS ###\n");
    if (use_DST) {
        delete dst;
    }
    else {
    }
    return res;
}
