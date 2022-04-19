#include <bits/stdc++.h>
#include "dst.h"
#include "util.h"

#define BLOOM_MAX_LAYER 63
#define CUCKOO_MAX_LAYER 63
#define DATASET_PATH std::string("/home/gt/Dataset/SOSD/data/")
#define CONFIG_PATH std::string("/home/gt/Dataset/SOSD/fpr_opt/")

FILE* file;
vector<uint64_t> key;
vector<uint16_t> runid, lcp;
vector<Bitwise> key_bit;
vector<pair<uint64_t, size_t>> query;

uint64_t randu(){
    return (uint64_t) rand() << 62 ^ (uint64_t) rand() << 31 ^ (uint64_t) rand();
}

void run(double totMem, string dataset_name, uint64_t cuckoo_mask){
    clock_t begin, end;
    size_t io = 0;
    function<pair<vector<size_t>, vector<size_t>> (vector<size_t>, vector<size_t>, size_t, size_t, uint64_t)> func;
    vector<double> bpn_bf, bpn_ck;
    {
        std::ifstream file((CONFIG_PATH + to_string((int)round(totMem / key.size())) + '_' + dataset_name).c_str(), std::ifstream::in);
        string token;
        for(bool bf = true; file >> token;)
            if(token == "@bpn"){
                file >> token >> token;
                (bf ? bpn_bf : bpn_ck).push_back(stof(token));
                bf = !bf;
            }
        file.close();
    }
    func = [totMem, &bpn_bf, &bpn_ck](vector<size_t> idist, vector<size_t> ldist, size_t bf_max, size_t ck_max, uint64_t ck_mask) -> pair<vector<size_t>, vector<size_t>> {
        for(int i = 0; i < 64; i++){
            idist[i] = max(idist[i] * bpn_bf[i + 1], 1.);
            ldist[i] = max(ldist[i] * bpn_ck[i + 1], 1.);
        }
        return make_pair(idist, ldist);
    };
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

    SplittedRosetta<vacuum::VacuumFilter<uint32_t>> dst(64, BLOOM_MAX_LAYER, CUCKOO_MAX_LAYER, cuckoo_mask, func, io_sim);

    fprintf(stderr, "Start Insert\n");
    begin = clock();
    if(!dst.AddKeys(key_bit, runid))
        return;
    end = clock();
    fprintf(stderr, "Finish Insert\n");
    fprintf(file, "InsertTP %lf  ", (double)key.size()*CLOCKS_PER_SEC/1e6/(end-begin));

    begin = clock();
    for (size_t i = 0; i < query.size(); i++) {
        from = query[i].first;
        ans = query[i].second;
        correct = false;
        auto res = dst.Seek(Bitwise(from));
        if(!correct){
            printf("WRONG ANSWER! query id: %lu;  from: %lu, ans: %lu\n", i, from, key[ans]);
            throw std::runtime_error("False Negative");
            return;
        }
        delete res;
    }
    end = clock();

    fprintf(file, "QueryTP %lf  ", (double)query.size()*CLOCKS_PER_SEC/1e6/(end-begin));
    fprintf(file, "IO %lf  ", (double)io/query.size());
    fprintf(file, "BPK %lf  ", (double)dst.mem() * 8 / key.size());
    fprintf(file, "expectedBPK %lf  ", (double)totMem / key.size());
    fprintf(file, "\n");
    fflush(file);
}
void test(string filename_prefix, string dataset_name, size_t n, size_t run_n, double bpkMin, double bpkMax, uint64_t cuckoo_mask){
    // Obtain a BPK - I/O cost image

    file = filename_prefix == "" ? stdout : fopen(("./log/" + filename_prefix + dataset_name + ".txt").c_str(), "w");

    size_t q = 1e6;
    key.clear();
    runid.clear();
    key_bit.clear();
    lcp.clear();
    query.clear();
#ifdef DATASET_PATH
    key = util::load_data<uint64_t>(DATASET_PATH + dataset_name);
#else
    for (size_t i=0; i<n; ++i)
        key.push_back(randu());
    sort(key.begin(), key.end());
#endif
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
        run(bpk * key.size(), dataset_name, cuckoo_mask);
        fprintf(stderr, "Done [n = %lu, bpk = %lf]\n", n, bpk);
        step = 2; // step = bpk * (pow(2, 0.1) - 1); // set step
    }
    if(file != stdout)
        fclose(file);
}

int main() {
    srand(1234);
    std::filesystem::create_directory("log");

    test(std::string("test-"), "books_200M_uint64", 2e8, 100, 18, 22, 0xfffffffff8000000ull);
    // test(std::string("4-18-"), "books_200M_uint64", 2e8, 100, 14, 22, 0xfffffffff8000000ull);
    // test(std::string("4-18-"), "normal_200M_uint64", 2e8, 100, 14, 22, 0xfffffffffc000000ull);
    // test(std::string("4-18-"), "uniform_sparse_200M_uint64", 2e8, 100, 14, 22, 0xffffffffff000000ull);
    // test(std::string("4-18-"), "osm_cellids_200M_uint64", 2e8, 100, 14, 22, 0xffffffff80000000ull);
    // test(std::string("4-18-"), "lognormal_200M_uint64", 2e8, 100, 14, 22, 0xffffffc000000000ull);

    return 0;
}
