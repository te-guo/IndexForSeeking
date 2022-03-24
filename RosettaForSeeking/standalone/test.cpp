#include <bits/stdc++.h>
#include "dst.h"

#define UNLAYERED
#define CUCKOO_FP_LEN 4

FILE* file;
vector<uint32_t> in_vec;
vector<Bitwise> in_bit;
vector<pair<uint32_t, uint32_t>> query;

void run(size_t diffidence, double memScale, int funcType, double funcPara){
    clock_t begin, end;
    size_t io = 0;
#ifdef UNLAYERED
    function<size_t (vector<size_t>)> func = [memScale](vector<size_t> distribution) -> size_t {
        size_t mem = 0;
        for(auto &i: distribution)
            mem += i;
        return ceil(mem * memScale);
    };
#else
    function<vector<size_t> (vector<size_t>)> func;
    if(funcType == 0)
        func = [memScale](vector<size_t> distribution) -> vector<size_t> {
            for(auto &i: distribution)
                i = ceil(i*memScale);
            return distribution;
        }; // bpn = const
    else if(funcType==1)
        func = [memScale, funcPara](vector<size_t> distribution) -> vector<size_t> {
            double k;
            double tot = 0, used = 0;
            k = 1;
            for(auto &i: distribution){
                k *= funcPara;
                tot += i;
                used += i * k;
            }
            k = tot * memScale / used;
            for(auto &i: distribution){
                k *= funcPara;
                i = ceil(i * k);
            }
            return distribution;
        }; // bpn = pow(para, -h)
    else if(funcType==2)
        func = [memScale, funcPara](vector<size_t> distribution) -> vector<size_t> {
            double tot = 0, k;
            for(auto &i: distribution)
                tot += i;
            k = tot * memScale / (tot + distribution.back() * funcPara);
            for(auto i = distribution.begin(); i != distribution.end(); i++)
                if(i + 1 == distribution.end())
                    *i = ceil(*i * k * funcPara);
                else
                    *i = ceil(*i * k);
            return distribution;
        }; // bpn = h == 0 ? para : 1
    else if(funcType==3)
        func = [memScale, funcPara](vector<size_t> distribution) -> vector<size_t> {
            double tot = 0, used = 0, k;
            for(auto i = distribution.begin(); i != distribution.end(); i++){
                tot += *i;
                used += *i * pow(funcPara, 1 / (distribution.end() - i));
            }
            k = tot * memScale / used;
            for(auto i = distribution.begin(); i != distribution.end(); i++)
                *i = ceil(*i * pow(funcPara, 1 / (distribution.end() - i)) * k);
            return distribution;
        }; // bpn = pow(para, 1 / (h + 1))
#endif

#ifdef UNLAYERED
    UnlayeredRosetta<
#else
    Rosetta<
#endif
#if CUCKOO_FP_LEN != 0
    vacuum::VacuumFilter<uint32_t, CUCKOO_FP_LEN + 1>
#else
    BloomFilter<>
#endif
    > dst(diffidence, 32, func);

    begin = clock();
    if(!dst.AddKeys(in_bit))
        return;
    end = clock();
    fprintf(file, "InsertTP %lf  ", (double)in_bit.size()*CLOCKS_PER_SEC/1e6/(end-begin));

    begin = clock();
    for (size_t i = 0; i < query.size(); i++) {
        uint32_t from = query[i].first, ans = query[i].second;
        while(true){
            uint32_t res = dst.Seek(Bitwise(from))->to_uint64()>>32;
            ++io;
            if (res == ans)
                break;
            else if (res < from || res > ans){
                printf("WRONG ANSWER! query id: %lu;  from: %u, res: %u, ans: %u\n", i, from, res, ans);
                throw std::runtime_error("False Negative");
                return;
            }
            from = res + 1;
        }
    }
    end = clock();

    //dst.printStats();
    fprintf(file, "QueryTP %lf  ", (double)query.size()*CLOCKS_PER_SEC/1e6/(end-begin));
    fprintf(file, "IO %lf  ", (double)io/query.size());
    fprintf(file, "BPK %lf  ", (double)dst.mem() * 8 / in_vec.size());
#ifdef UNLAYERED
    fprintf(file, "LF %lf  ", dst.get_load_factor());
#endif
    fprintf(file, "\n");
}
void test(string filename, size_t n, size_t diffidence, double scaleMin, double scaleMax, int funcType = 0, double funcPara = 0){
    // Obtain a BPK - I/O cost image

    file = filename == "" ? stdout : fopen(("./log/" + filename + ".txt").c_str(), "w");

    size_t q = 1e6;
    in_vec.clear();
    in_bit.clear();
    query.clear();
    for (size_t i=0; i<n; ++i)
        in_vec.push_back(rand());
    sort(in_vec.begin(), in_vec.end());
    for (size_t i=0; i<in_vec.size(); ++i)
        in_bit.push_back(in_vec[i]);
    for (size_t i=0; i<q; ++i) {
        uint32_t from = rand();
        auto it = lower_bound(in_vec.begin(), in_vec.end(), from);
        if(it == in_vec.end()){
            --i;
            continue;
        }
        query.push_back(make_pair(from, *it));
    }

    fprintf(file, "Insert Throughput (M/s), Query Throughput (M/s), Expected I/O Cost, BPK\nEvaluation:\n");
    for(double scale=scaleMin, step=0.1; scale <= scaleMax; scale += step){
        run(diffidence, scale, funcType, funcPara);
        fprintf(stderr, "%lu %lu %lf %lf\n", n, diffidence, funcPara, scale);
        step = scale * (pow(2, 0.1) - 1); // set step
    }
    if(file != stdout)
        fclose(file);
}

int main() {
    srand(1234);
    std::filesystem::create_directory("log");

    // test("bloom_1e6", 1e6, 100, 1.5, 6, 2, 7);
    // test("bloom_1e7", 1e7, 100, 1.5, 6, 2, 7);
    // test("cuckoo_4bit_1e6", 1e6, 100, 4, 20, 0, 0);
    // test("cuckoo_8bit_1e6", 1e6, 100, 8, 20, 0, 0);
    // test("cuckoo_16bit_1e6", 1e6, 100, 16, 20, 0, 0);
    // test("unlayered_4bit_1e6", 1e6, 100, 4, 20, 0, 0);
    // test("unlayered_8bit_1e6", 1e6, 100, 8, 40, 0, 0);
    // test("unlayered_16bit_1e6", 1e6, 100, 16, 80, 0, 0);
    test("unlayered_4bit_5e7", 5e7, 100, 4, 20, 0, 0);
    // test("unlayered_8bit_5e7", 5e7, 100, 8, 40, 0, 0);
    // test("unlayered_16bit_5e7", 5e7, 100, 16, 80, 0, 0);

    return 0;
}
