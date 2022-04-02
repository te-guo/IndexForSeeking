#include <bits/stdc++.h>
#include "dst.h"

#define CUCKOO_FP_LEN 4
#define CUCKOO_MASK 0xfffffffff0000000ull

FILE* file;
vector<uint64_t> in_vec;
vector<Bitwise> in_bit;
vector<pair<uint64_t, pair<uint64_t, int>>> query;

uint64_t randu(){
    return (uint64_t) rand() << 62 ^ (uint64_t) rand() << 31 ^ (uint64_t) rand();
}

void run(double totMem, int funcType, double funcPara){
    clock_t begin, end;
    size_t io = 0;
    function<pair<vector<size_t>, vector<size_t>> (vector<size_t>, vector<size_t>, size_t, size_t, uint64_t)> func = [totMem](vector<size_t> idist, vector<size_t> ldist, size_t bf_max, size_t ck_max, uint64_t ck_mask) -> pair<vector<size_t>, vector<size_t>> {
        for(auto &i: ldist) i *= 3;
        for(int i=bf_max+1; i<64; i++)
            if(i > bf_max){
                idist[bf_max] += idist[i];
                idist[i] = 0;
            }
        for(int i=ck_max+1; i<64; i++)
            if(i > ck_max){
                ldist[ck_max] += ldist[i];
                ldist[i] = 0;
            }
        size_t mem = 0;
        for(auto &i: idist) mem += i;
        for(auto &i: ldist) mem += i;
        for(auto &i: idist) i = ceil(i * totMem / mem);
        for(auto &i: ldist) i = ceil(i * totMem / mem);
        return make_pair(idist, ldist);
    };
    size_t ans, ans_len;
    bool correct;
    function<bool (const Bitwise&)> io_sim = [&ans, &ans_len, &io, &correct](const Bitwise& str) -> bool {
        ++io;
        return correct |= ans_len == str.size() && ans >> 64 - ans_len == str.to_uint64() >> 64 - ans_len;
    };

    SplittedRosetta<vacuum::VacuumFilter<uint32_t, CUCKOO_FP_LEN + 1>> dst(64, 35, 28, CUCKOO_MASK, func, io_sim);

    begin = clock();
    if(!dst.AddKeys(in_bit))
        return;
    end = clock();
    std::cerr<<"Finish Insert\n";
    fprintf(file, "InsertTP %lf  ", (double)in_bit.size()*CLOCKS_PER_SEC/1e6/(end-begin));

    begin = clock();
    for (size_t i = 0; i < query.size(); i++) {
        uint64_t from = query[i].first;
        ans = query[i].second.first;
        ans_len = query[i].second.second;
        correct = false;
        auto res = dst.Seek(Bitwise(from));
        if(!correct){
            printf("WRONG ANSWER! query id: %lu;  from: %lu, ans: %lu\n", i, from, ans);
            throw std::runtime_error("False Negative");
            return;
        }
        delete res;
    }
    end = clock();

    fprintf(file, "QueryTP %lf  ", (double)query.size()*CLOCKS_PER_SEC/1e6/(end-begin));
    fprintf(file, "IO %lf  ", (double)io/query.size());
    fprintf(file, "BPK %lf  ", (double)dst.mem() * 8 / in_vec.size());
    fprintf(file, "\n");
    fflush(file);
}
void test(string filename, size_t n, double bpkMin, double bpkMax, int funcType = 0, double funcPara = 0){
    // Obtain a BPK - I/O cost image

    file = filename == "" ? stdout : fopen(("./log/" + filename + ".txt").c_str(), "w");

    size_t q = 1e6;
    in_vec.clear();
    in_bit.clear();
    query.clear();
    for (size_t i=0; i<n; ++i)
        in_vec.push_back(randu());
    sort(in_vec.begin(), in_vec.end());
    in_vec.erase(unique(in_vec.begin(), in_vec.end()), in_vec.end());
    for (size_t i=0; i<in_vec.size(); ++i)
        in_bit.push_back(in_vec[i]);
    for (size_t i=0; i<q; ++i) {
        uint64_t from = randu();
        auto it = lower_bound(in_vec.begin(), in_vec.end(), from) - in_vec.begin();
        if(it == in_vec.size()){
            --i;
            continue;
        }
        int unique_len = 0;
        if(it > 0)
            unique_len = max(unique_len, __builtin_clzll(in_vec[it] ^ in_vec[it-1]));
        if(it < in_vec.size() - 1)
            unique_len = max(unique_len, __builtin_clzll(in_vec[it] ^ in_vec[it+1]));
        while(!(CUCKOO_MASK >> unique_len & 1ull))
            unique_len++;
        query.push_back(make_pair(from, make_pair(in_vec[it], unique_len + 1)));
    }

    fprintf(file, "Insert Throughput (M/s), Query Throughput (M/s), Expected I/O Cost, BPK\nEvaluation:\n");
    for(double bpk=bpkMin, step=1; bpk <= bpkMax + 1e-6; bpk += step){
        run(bpk * in_vec.size(), funcType, funcPara);
        fprintf(stderr, "%lu %lf %lf\n", n, funcPara, bpk);
        // step = bpk * (pow(2, 0.1) - 1); // set step
    }
    if(file != stdout)
        fclose(file);
}

int main() {
    srand(1234);
    std::filesystem::create_directory("log");

    test("", 2e7 /*1e8*/, 20, 40, 0, 0);

    return 0;
}
