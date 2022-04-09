#include <bits/stdc++.h>
#include "dst.h"

#define CUCKOO_MASK 0xfffffffff0000000ull
#define CUCKOO_MAX_LAYER 28
#define BLOOM_MAX_LAYER 34
#define GUROBI
// #define UNSPLITTED

FILE* file;
vector<uint64_t> key;
vector<uint16_t> runid, lcp;
vector<Bitwise> key_bit;
vector<pair<uint64_t, size_t>> query;

uint64_t randu(){
    return (uint64_t) rand() << 62 ^ (uint64_t) rand() << 31 ^ (uint64_t) rand();
}

void run(double totMem, int funcType, double funcPara){
    clock_t begin, end;
    size_t io = 0;
#ifndef UNSPLITTED
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
#ifdef GUROBI
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
    int idx = round((totMem / key.size() - 14) / 2);
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
    uint64_t from;
    size_t ans;
    bool correct;
    function<bool (uint16_t, const Bitwise&)> io_sim = [&](uint16_t rid, const Bitwise& prefix) -> bool {
        ++io;
        // should be: find the only key in the rid-th run that has the expected prefix, and check the key >= from.
        if(rid != runid[ans]) return false;
        if(prefix.lcp(Bitwise(key_bit[ans], prefix.size())) < prefix.size()) return false;
        if(prefix.size() <= lcp[ans]) return false;
        return correct = true;
    };

    SplittedRosetta<vacuum::VacuumFilter<uint32_t>> dst(64, BLOOM_MAX_LAYER, CUCKOO_MAX_LAYER, CUCKOO_MASK, func, io_sim);
#else
    function<vector<size_t> (vector<size_t>)> func;
    if(funcType == 0)
        func = [totMem](vector<size_t> distribution) -> vector<size_t> {
            size_t mem = 0;
            for(auto &i: distribution) mem += i;
            for(auto &i: distribution) i = ceil(i * totMem / mem);
            return distribution;
        }; // bpn = const
    else if(funcType==1)
        func = [totMem, funcPara](vector<size_t> distribution) -> vector<size_t> {
            double k = 1;
            for(auto &i: distribution){
                k *= funcPara;
                i = ceil(i * k);
            }
            size_t mem = 0;
            for(auto &i: distribution) mem += i;
            for(auto &i: distribution) i = ceil(i * totMem / mem);
            return distribution;
        }; // bpn = pow(para, -h)
    else if(funcType==2)
        func = [totMem, funcPara](vector<size_t> distribution) -> vector<size_t> {
            distribution.back() = ceil(distribution.back() * funcPara);
            size_t mem = 0;
            for(auto &i: distribution) mem += i;
            for(auto &i: distribution) i = ceil(i * totMem / mem);
            return distribution;
        }; // bpn = h == 0 ? para : 1

    Rosetta<
#ifdef CUCKOO_FP_LEN
    vacuum::VacuumFilter<uint32_t, CUCKOO_FP_LEN + 1>
#else
    BloomFilter<>
#endif
    > dst(500, 25, func);
#endif

    fprintf(stderr, "Start Insert\n");
    begin = clock();
    if(!dst.AddKeys(key_bit, runid))
        return;
    end = clock();
    fprintf(stderr, "Finish Insert\n");
    fprintf(file, "InsertTP %lf  ", (double)key.size()*CLOCKS_PER_SEC/1e6/(end-begin));

    begin = clock();
#ifndef UNSPLITTED
    for (size_t i = 0; i < query.size(); i++) {
        from = query[i].first;
        ans = query[i].second;
        correct = false;
        auto res = dst.Seek(Bitwise(from));
        if(!correct){
            printf("WRONG ANSWER! query id: %lu;  from: %lu, ans: %lu\n", i, from, ans);
            throw std::runtime_error("False Negative");
            return;
        }
        delete res;
    }
#else
    for (size_t i = 0; i < query.size(); i++) {
        uint64_t from = query[i].first, ans = query[i].second.first;
        while(true){
            auto out = dst.Seek(Bitwise(from));
            uint64_t res = out->to_uint64();
            ++io;
            if (res == ans)
                break;
            else if (res < from || res > ans){
                printf("WRONG ANSWER! query id: %lu;  from: %lu, res: %lu, ans: %lu\n", i, from, res, ans);
                throw std::runtime_error("False Negative");
                return;
            }
            from = res + 1;
            delete out;
        }
    }
#endif
    end = clock();

    fprintf(file, "QueryTP %lf  ", (double)query.size()*CLOCKS_PER_SEC/1e6/(end-begin));
    fprintf(file, "IO %lf  ", (double)io/query.size());
    fprintf(file, "BPK %lf  ", (double)dst.mem() * 8 / key.size());
    fprintf(file, "\n");
    fflush(file);
}
void test(string filename, size_t n, size_t run_n, double bpkMin, double bpkMax, int funcType = 0, double funcPara = 0){
    // Obtain a BPK - I/O cost image

    file = filename == "" ? stdout : fopen(("./log/" + filename + ".txt").c_str(), "w");

    size_t q = 1e6;
    key.clear();
    runid.clear();
    key_bit.clear();
    lcp.clear();
    query.clear();
    for (size_t i=0; i<n; ++i)
        key.push_back(randu());
    sort(key.begin(), key.end());
    key.erase(unique(key.begin(), key.end()), key.end());
    for (size_t i=0; i<key.size(); ++i)
        runid.push_back(rand() % run_n);
    for (size_t i=0; i<key.size(); ++i)
        key_bit.push_back(key[i]);
    for (size_t i=0, l1=0; i<key.size(); ++i){
        size_t l2 = i + 1 < key.size() ? key_bit[i].lcp(key_bit[i + 1]) : 0;
        lcp.push_back(max(l1, l2));
        l1 = l2;
    }
    for (size_t i=0; i<q; ++i) {
        uint64_t from = randu();
        auto it = lower_bound(key.begin(), key.end(), from) - key.begin();
        if(it == key.size()){
            --i;
            continue;
        }
        query.push_back(make_pair(from, it));
    }

    fprintf(file, "Insert Throughput (M/s), Query Throughput (M/s), Expected I/O Cost, BPK\nEvaluation:\n");
    for(double bpk=bpkMin, step=1; bpk <= bpkMax + 1e-6; bpk += step){
        run(bpk * key.size(), funcType, funcPara);
        fprintf(stderr, "Done [n = %lu, funcPara = %lf, bpk = %lf]\n", n, funcPara, bpk);
        step = 2; // step = bpk * (pow(2, 0.1) - 1); // set step
    }
    if(file != stdout)
        fclose(file);
}

int main() {
    srand(1234);
    std::filesystem::create_directory("log");

    test(""/*"splitted_5e8_l34"*/, 5e8, 100, 14, 22, 0, 0);

    return 0;
}
