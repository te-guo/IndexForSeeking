#include <bits/stdc++.h>
#include "dst.h"

#define CUCKOO_FP_LEN 4
#define CUCKOO_MASK 0xfffffffff0000000ull
#define CUCKOO_MAX_LAYER 28
#define BLOOM_MAX_LAYER 34
#define gurobi

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
#ifdef gurobi
    static double gurobi_a[][7]={
        {0.3551073730976299, 1.599990748245976, 0.6152689257681739, 0.6072477971431384, 0.4176961501265128, 0.42144577802115324, 0.508743676268198},
        {0.29465396937148214, 0.09208012199783208, 1.6, 0.9808025571881682, 0.1951027573978623, 0.44560027309238026, 0.322239115845548},
        {0.04843458434602678, 1.0192572784259177, 1.0008917302035203, 0.09417096671082589, 0.13169989844310842, 0.05383147342764956, 1.0008917302035205},
        {0.38538357672813534, 0.1575087242395704, 0.06704623108234742, 1.5999999999999828, 0.9428213939112957, 0.2995820480390141, 0.06689184990075958},
        {0.10989799699907793, 0.05540779799357534, 1.4644447693046942, 0.05383147342764956, 0.12043013123591745, 0.10124884290074966, 0.08419979891931928}
    };
    static double gurobi_b[]={
        0.002350783115917231,
        0.0010024389546877874,
        0.0009538004771047486,
        0.0004492774983996558,
        0.0002266338277448154
    };
    int idx = round((totMem / in_vec.size() - 14) / 2);
    func = [totMem, idx](vector<size_t> idist, vector<size_t> ldist, size_t bf_max, size_t ck_max, uint64_t ck_mask) -> pair<vector<size_t>, vector<size_t>> {
        for(int i=bf_max+1; i<64; i++)
            if(i > bf_max){
                idist[bf_max] += idist[i];
                idist[i] = 0;
            }
        for(int i=bf_max, j=6; j>=0; i--, j--)
            idist[i] = ceil(idist[i] * (-log(gurobi_a[idx][j]/2)/pow(log(2),2)));
        for(int i=bf_max-7, j=1; i>=0; i--, j++)
            idist[i] = ceil(idist[i] * (-log(gurobi_a[idx][0]/2)/pow(log(2),2)/pow(1.25, j)));
        for(int i=ck_max+1; i<64; i++)
            if(i > ck_max){
                ldist[ck_max] += ldist[i];
                ldist[i] = 0;
            }
        ldist[ck_max] = ceil(ldist[ck_max] * (log(4/gurobi_b[idx])/log(2)));
        size_t mem = 0;
        for(auto &i: idist) mem += i;
        for(auto &i: ldist) mem += i;
        for(auto &i: idist) i = ceil(i * totMem / mem);
        for(auto &i: ldist) i = ceil(i * totMem / mem);
        return make_pair(idist, ldist);
    };
#endif
    size_t ans, ans_len;
    bool correct;
    function<bool (const Bitwise&)> io_sim = [&ans, &ans_len, &io, &correct](const Bitwise& str) -> bool {
        ++io;
        return correct |= ans_len == str.size() && ans >> 64 - ans_len == str.to_uint64() >> 64 - ans_len;
    };

    SplittedRosetta<vacuum::VacuumFilter<uint32_t, CUCKOO_FP_LEN + 1>> dst(64, BLOOM_MAX_LAYER, CUCKOO_MAX_LAYER, CUCKOO_MASK, func, io_sim);

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
        step = 2; // step = bpk * (pow(2, 0.1) - 1); // set step
    }
    if(file != stdout)
        fclose(file);
}

int main() {
    srand(1234);
    std::filesystem::create_directory("log");

    test("splitted_4bit_5e8", 5e8, 14, 22, 0, 0);

    return 0;
}
