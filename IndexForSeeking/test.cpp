#include <bits/stdc++.h>
#include "dst.h"
#include "util.h"

#define DATASET_PATH std::string("/home/gt/Dataset/SOSD/data/")
#define CONFIG_PATH std::string("/home/gt/Dataset/SOSD/fpr_offset/")
// #define ZIPFIAN_DIST
#ifdef ZIPFIAN_DIST
#define ZIPF_ALPHA 1.5
#define ZIPF_GROUP 1000000
#endif

FILE* file;
vector<uint64_t> key;
vector<uint16_t> runid, lcp;
vector<Bitwise> key_bit;
vector<pair<uint64_t, size_t>> query;

uint64_t randu(){
    return (uint64_t) rand() << 62 ^ (uint64_t) rand() << 31 ^ (uint64_t) rand();
}

void run(double totMem, string dataset_name){
    clock_t begin, end;
    size_t io = 0, io_hash = 0;
#ifdef ZIPFIAN_DIST
    vector<size_t> ios;
    size_t last_io = 0;
#endif
    function<pair<vector<size_t>, vector<size_t>> (vector<size_t>, vector<size_t>)> func;
    vector<double> bpn_bf, bpn_ck;
    size_t shift = 0, padding = 0;
    {
        std::ifstream file((CONFIG_PATH + to_string((int)round(totMem / key.size())) + '_' + dataset_name).c_str(), std::ifstream::in);
        string token;
        for(bool bf = true; file >> token;)
            if(token == "@bpn"){
                file >> token; if(token != "=") continue;
                file >> token;
                (bf ? bpn_bf : bpn_ck).push_back(stof(token));
                bf = !bf;
            }
            else if(token == "Shifting"){
                file >> token; if(token != "Down:") continue;
                file >> token;
                shift = atoi(token.c_str());
            }
            else if(token == "Padding:"){
                file >> token;
                padding = atoi(token.c_str()) - 1;
            }
        file.close();
    }
    func = [totMem, &bpn_bf, &bpn_ck](vector<size_t> idist, vector<size_t> ldist) -> pair<vector<size_t>, vector<size_t>> {
        for(int i = 0; i < 64; i++){
            idist[i] = ceil(idist[i] * bpn_bf[i + 1]);
            ldist[i] = ceil(ldist[i] * bpn_ck[i + 1]);
        }
        return make_pair(idist, ldist);
    };
    uint64_t from;
    size_t ans;
    bool correct;
    function<int (uint16_t, const Bitwise&)> io_sim = [&](uint16_t rid, const Bitwise& prefix) -> int {
        ++io;
        // should be: find the only key in the rid-th run that has the expected prefix, and check the key >= from.
        if(rid == runid[ans]
            && prefix.lcp(Bitwise(key_bit[ans], prefix.size())) >= prefix.size()
            && prefix.size() > lcp[ans]) return correct = 1;
        if(ans > 0 && rid == runid[ans - 1]
            && prefix.lcp(Bitwise(key_bit[ans - 1], prefix.size())) >= prefix.size()
            && prefix.size() > lcp[ans - 1]) return -1;
        ++io_hash;
        return 0;
    };

    SplittedRosetta<vacuum::VacuumFilter<uint32_t>> dst(64, shift, padding, func, io_sim);

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
#ifdef ZIPFIAN_DIST
        if(i % ZIPF_GROUP == ZIPF_GROUP - 1){
            ios.push_back(io - last_io);
            last_io = io;
        }
#endif
    }
    end = clock();

#ifdef ZIPFIAN_DIST
    {
        double avg = 0, sd = 0;
        for(auto &x: ios)
            avg += (double)x / ZIPF_GROUP;
        avg /= ios.size();
        for(auto &x: ios)
            sd += ((double)x / ZIPF_GROUP - avg) * ((double)x / ZIPF_GROUP - avg);
        sd = sqrt(sd / ios.size());
        fprintf(file, "SD_IO %lf  ", sd);
    }
#endif
    fprintf(file, "QueryTP %lf  ", (double)query.size()*CLOCKS_PER_SEC/1e6/(end-begin));
    fprintf(file, "IO %lf  ", (double)io/query.size());
    fprintf(file, "IO_hash %lf  ", (double)io_hash/query.size());
    fprintf(file, "BPK %lf  ", (double)dst.mem() * 8 / key.size());
    fprintf(file, "expectedBPK %lf  ", (double)totMem / key.size());
    fprintf(file, "\n");
    fflush(file);
}
void test(string filename_prefix, string dataset_name, size_t n, size_t run_n, double bpkMin, double bpkMax){
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
#ifdef ZIPFIAN_DIST
    zipf::init(ZIPF_ALPHA, q);
    map<int, pair<uint64_t, size_t>> zipf_map;
    q = q / ZIPF_GROUP * ZIPF_GROUP;
#endif
    for (size_t i=0; i<q; ++i) {
        uint64_t from;
        size_t it;
#ifndef ZIPFIAN_DIST
        from = key.back() == 0xffffffffffffffffULL ? randu() : randu() % (key.back() + 1ULL);
        it = lower_bound(key.begin(), key.end(), from) - key.begin();
#else
        int idx = zipf::zipf();
        if(!zipf_map.count(idx)){
            from = key.back() == 0xffffffffffffffffULL ? randu() : randu() % (key.back() + 1ULL);
            it = lower_bound(key.begin(), key.end(), from) - key.begin();
            zipf_map[idx] = make_pair(from, it);
        }
        else
            from = zipf_map[idx].first, it = zipf_map[idx].second;
#endif
        query.push_back(make_pair(from, it));
#ifdef ZIPFIAN_DIST
        if(i % ZIPF_GROUP == ZIPF_GROUP - 1)
            zipf_map.clear();
#endif
    }

    fprintf(file, "Insert Throughput (M/s), Query Throughput (M/s), Expected I/O Cost, BPK\nEvaluation:\n");
    for(double bpk=bpkMin, step=1; bpk <= bpkMax + 1e-6; bpk += step){
        run(bpk * key.size(), dataset_name);
        fprintf(stderr, "Done [n = %lu, bpk = %lf]\n", n, bpk);
        step = 2; // step = bpk * (pow(2, 0.1) - 1); // set step
    }
    if(file != stdout)
        fclose(file);
}

int main() {
    srand(1234);
    std::filesystem::create_directory("log");

    // test(std::string("4-21-"), "books_200M_uint64", 2e8, 100, 12, 22);
    // test(std::string("4-21-"), "normal_200M_uint64", 2e8, 100, 12, 22);
    // test(std::string("4-21-"), "uniform_sparse_200M_uint64", 2e8, 100, 12, 22);
    // test(std::string("4-21-"), "osm_cellids_200M_uint64", 2e8, 100, 12, 22);
    // test(std::string("4-21-"), "lognormal_200M_uint64", 2e8, 100, 12, 22);
    // test(std::string("4-21-"), "fb_200M_uint64", 2e8, 100, 12, 22);
    // test(std::string("4-21-"), "wiki_ts_200M_uint64", 2e8, 100, 12, 22);
    // test(std::string("4-21-"), "uniform_dense_200M_uint64", 2e8, 100, 12, 22);

    return 0;
}
