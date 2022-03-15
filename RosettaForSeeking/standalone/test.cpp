#include <bits/stdc++.h>
#include "dst.h"

FILE* file;
vector<uint32_t> in_vec;
vector<Bitwise> in_bit;
vector<pair<uint32_t, uint32_t>> query;

void run(size_t diffidence, double memScale, int funcType, double funcPara){
    clock_t begin, end;
    size_t io = 0;
    function<vector<size_t> (vector<size_t>)> func;

    if(funcType == 0)
        func = [memScale](vector<size_t> distribution) -> vector<size_t> {
            for(auto &i: distribution)
                i = ceil(i*memScale);
            return distribution;
        };
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
        };
    DstFilter<BloomFilter<>> dst(diffidence, 32, func);

    begin = clock();
    dst.AddKeys(in_bit);
    end = clock();
    fprintf(file, "InsertTP %lf  ", (double)in_bit.size()*CLOCKS_PER_SEC/1e6/(end-begin));

    begin = clock();
    for (auto &q: query) {
        uint32_t from = q.first, ans = q.second;
        while(true){
            uint32_t res = ntohl(*(uint32_t*)(dst.Seek(Bitwise(from))->data()));
            ++io;
            if (res == ans)
                break;
            else if (res < from || res > ans){
                printf("WRONG ANSWER!  from: %u, res: %u, ans: %u\n", from, res, ans);
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
        if(scale<3)step=0.1; else if(scale<5)step=0.15; else if(scale<7)step=0.25; else step=0.4; // set step
    }
    if(file != stdout)
        fclose(file);
}

int main() {
    std::filesystem::create_directory("log");

/*
    test("d=50", 1e5, 50, 1, 2, 0, 0);
    test("d=100", 1e5, 100, 1, 2, 0, 0);
    test("d=300", 1e5, 300, 1, 2, 0, 0);
    test("d=1000", 1e5, 1000, 1, 2, 0, 0);
    test("d=3000", 1e5, 3000, 1, 2, 0, 0);
*/

/*
    test("f=const", 1e6, 100, 2, 5, 0, 0);
    test("f=pow_1.25", 1e6, 100, 2, 5, 1, 1.25);
    test("f=pow_1.35", 1e6, 100, 2, 5, 1, 1.35);
    test("f=pow_1.15", 1e6, 100, 2, 5, 1, 1.15);
*/

    test("N=1e4", 1e4, 100, 1.5, 5, 1, 1.25);
    test("N=1e5", 1e5, 100, 1.5, 5, 1, 1.25);
    test("N=1e6", 1e6, 100, 1.5, 5, 1, 1.25);
    test("N=1e7", 1e7, 100, 1.5, 5, 1, 1.25);

    return 0;
}
