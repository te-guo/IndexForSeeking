#ifndef DST_H_
#define DST_H_

#include <iostream>
#include <algorithm>
#include <vector>
#include <memory>
#include <cassert>
#include <functional>
#include <arpa/inet.h>
#include <cstring>
#include <string>
#include "MurmurHash3.h"

#ifdef USE_DTL
#include <dtl/filter/blocked_bloomfilter/blocked_bloomfilter.hpp>
#include <dtl/filter/blocked_bloomfilter/blocked_bloomfilter_tune.hpp>
#include <dtl/filter/bbf_64.hpp>
#endif

using namespace std;

vector<size_t> calc_dst(vector<size_t> dist, double bpk, vector<size_t> qdist, size_t cutoff = 0);

template<typename T>
class FastModulo {
    private:
        T num_;

    public:
        FastModulo(T num): num_(num) { }
        friend T operator%(const T &v, const FastModulo &u) {
            return v%u.num_;
        }
        T value() {
            return num_;
        }
};

class Bitwise {
    private:
        uint8_t *data_ = 0;
        size_t len_;
        bool is_mut_; // Create immutable views when copying data that will
                      // last less than the mutable version. Much faster than shared_ptr. A
                      // smart person would have done this through inheritance, but I am not
                      // a smart person.

    public:
        Bitwise(const string &key): len_(8*(key.size()+1)), is_mut_(true) {
            if (len_ == 0) return;
            data_ = new uint8_t[len_/8];
            assert(data_);
            memcpy(data_, key.c_str(), len_/8);
        }

        Bitwise(const uint32_t &key): len_(8*sizeof(uint32_t)), is_mut_(true) {
            if (len_ == 0) return;
            data_ = new uint8_t[len_/8];
            assert(data_);
            uint32_t tmp = htonl(key);
            memcpy(data_, &tmp, len_/8);
        }

        Bitwise(const uint64_t &key): len_(8*sizeof(uint64_t)), is_mut_(true) {
            if (len_ == 0) return;
            data_ = new uint8_t[len_/8];
            assert(data_);
            uint64_t tmp = __builtin_bswap64(key);
            memcpy(data_, &tmp, len_/8);
        }

        Bitwise(bool init, size_t len): len_(len), is_mut_(true) {
            if (len_ == 0) return;
            data_ = new uint8_t[(len_+7)/8]();
            assert(data_);
            if (init == 1) {
                assert(false);
                memset(data_, 255, (len_+7)/8);
            }
        }

        Bitwise(const uint8_t *data, size_t nbytes): len_(8*nbytes), is_mut_(true) {
            if (len_ == 0) return;
            data_ = new uint8_t[len_/8];
            assert(data_);
            memcpy(data_, data, len_/8);
        }

        Bitwise(uint8_t *data, size_t nbytes, bool is_mut): data_(data), len_(8*nbytes), is_mut_(is_mut) { }

        // Copy constructor should only be used when adding new levels to a non-empty DST
        Bitwise(const Bitwise &other): len_(other.size()), is_mut_(other.is_mut_) {
            assert(false);
            if (is_mut_) {
                data_ = new uint8_t[len_/8];
                assert(data_);
                memcpy(data_, other.data(), len_/8);
            }
            else
                data_ = other.data();
        }

        Bitwise(Bitwise &&other) noexcept {
            len_ = other.size();
            data_ = other.data();
            if (other.is_mut_) {
                other.is_mut_ = false;
                is_mut_ = true;
            }
            else
                is_mut_ = false;
        }

        // Needs to exist for some reason, but should never be called.
        Bitwise &operator=(const Bitwise &other) {
            assert(false);
            len_ = other.size();
            data_ = other.data();
            is_mut_ = false;
            return *this;
        }

        Bitwise(const Bitwise &other, size_t len): len_(len) {
            if (len <= other.size()) {
                data_ = other.data();
                is_mut_ = false;
            }
            else {
                data_ = new uint8_t[(len_+7)/8]();
                memcpy(data_, other.data(), other.size()/8);
                if (other.size()%8 != 0)
                    assert(false);
                is_mut_ = true;
            }
        }

        ~Bitwise() {
            if (is_mut_ and data_)
                delete[] data_;
        }

        uint64_t hash(uint32_t seed, FastModulo<uint64_t> modulo) const {
            if (modulo.value() < (1ULL<<32)) {
                uint32_t h;
                if (size() & 7u){
                    uint8_t tmp[size() + 7 >> 3];
                    copy_n(data(), size() + 7 >> 3, tmp);
                    tmp[size() - 1 >> 3] >>= 8 - (size() & 7u);
                    MurmurHash3_x86_32(tmp, size() + 7 >> 3, seed, &h);
                }
                else
                    MurmurHash3_x86_32(data(), size() >> 3, seed, &h);
                return h % modulo;
            }
            uint64_t h[2];
            if (size() & 7u){
                uint8_t tmp[size() + 7 >> 3];
                copy_n(data(), size() + 7 >> 3, tmp);
                tmp[size() - 1 >> 3] >>= 8 - (size() & 7u);
                MurmurHash3_x64_128(tmp, size() + 7 >> 3, seed, h);
            }
            else
                MurmurHash3_x64_128(data(), size() >> 3, seed, h);
            return h[0] % modulo;
        }
        uint64_t hash(uint32_t seed) const {
            uint64_t h[2];
            if (size() & 7u){
                uint8_t tmp[size() + 7 >> 3];
                copy_n(data(), size() + 7 >> 3, tmp);
                tmp[size() - 1 >> 3] >>= 8 - (size() & 7u);
                MurmurHash3_x64_128(tmp, size() + 7 >> 3, seed, h);
            }
            else
                MurmurHash3_x64_128(data(), size() >> 3, seed, h);
            return h[0];
        }
        uint64_t hash_len(uint32_t seed, FastModulo<uint64_t> modulo) const {
            if (modulo.value() < (1ULL<<32)) {
                uint32_t h;
                uint8_t tmp[size() + 15 >> 3];
                copy_n(data(), size() + 7 >> 3, tmp);
                tmp[size() >> 3] >>= 8 - (size() & 7u);
                tmp[size() + 7 >> 3] = len_;
                MurmurHash3_x86_32(tmp, size() + 15 >> 3, seed, &h);
                return h % modulo;
            }
            uint64_t h[2];
            uint8_t tmp[size() + 15 >> 3];
            copy_n(data(), size() + 7 >> 3, tmp);
            tmp[size() >> 3] >>= 8 - (size() & 7u);
            tmp[size() + 7 >> 3] = len_;
            MurmurHash3_x64_128(tmp, size() + 15 >> 3, seed, h);
            return h[0] % modulo;
        }
        uint64_t hash_len(uint32_t seed) const {
            uint64_t h[2];
            uint8_t tmp[size() + 15 >> 3];
            copy_n(data(), size() + 7 >> 3, tmp);
            tmp[size() >> 3] >>= 8 - (size() & 7u);
            tmp[size() + 7 >> 3] = len_;
            MurmurHash3_x64_128(tmp, size() + 15 >> 3, seed, h);
            return h[0];
        }

        uint64_t to_uint64() const {
            if (len_<=64) {
                uint64_t out = 0;
                memcpy(&out, data_, len_ + 7 >> 3);
                out = __builtin_bswap64(out);
                out >>= 64 - len_;
                out <<= 64 - len_;
                return out;
            }
            return 0; // TODO
        }

        void print() const {
            for (size_t i=0; i<len_; ++i)
                printf("%d", get(i));
        }

        bool get(size_t i) const {
            return (data_[i/8] >> (7-(i%8))) & 1;
        }

        void flip(size_t i) {
            assert(is_mut_);
            data_[i/8] ^= (1<<(7-(i%8)));
        }

        void set(size_t i, bool v) {
            assert(is_mut_);
            if (get(i) != v)
                flip(i);
        }

        size_t lcp(const Bitwise &other) const {
            size_t out = 0;
            while (out+64 <= other.size() && out+64 <= size() && ((uint64_t*)data_)[out/64] == ((uint64_t*)other.data())[out/64])
                out += 64;
            while (out+32 <= other.size() && out+32 <= size() && ((uint32_t*)data_)[out/32] == ((uint32_t*)other.data())[out/32])
                out += 32;
            while (out+16 <= other.size() && out+16 <= size() && ((uint16_t*)data_)[out/16] == ((uint16_t*)other.data())[out/16])
                out += 16;
            while (out+8 <= other.size() && out+8 <= size() && data_[out/8] == other.data()[out/8])
                out += 8;
            while (out+1 <= other.size() && out+1 <= size() && get(out) == other.get(out))
                out += 1;
            return out;
        }

        size_t size() const {
            return len_;
        }

        uint8_t *data() const {
            return data_;
        }
};

class Filter {
    public:
        Filter(){};
        virtual ~Filter(){};
        virtual bool AddKeys(const vector<Bitwise> &keys, const vector<uint16_t> & runids) = 0;
        virtual bool Query(const Bitwise &key) = 0;
};

template<bool keep_stats=false>
class BloomFilter final : public Filter {
    private:
        Bitwise data_;
        size_t nhf_;
        size_t nkeys_;
        vector<size_t> seeds_;
        FastModulo<uint64_t> nmod_;

    public:
        size_t nqueries_ = 0;
        size_t npositives_ = 0;

        BloomFilter(size_t nbits): data_(Bitwise(false, ((nbits+7)/8)*8)), nhf_(0), nkeys_(0), nmod_(((nbits+7)/8*8)) {}
        BloomFilter(size_t nbits, size_t nkeys): data_(Bitwise(false, ((nbits+7)/8)*8)), nhf_(0), nkeys_(nkeys), nmod_(((nbits+7)/8*8)) { init(); }
        BloomFilter(size_t nbits, uint8_t* data, size_t nhf, size_t nkeys, vector<size_t> seeds): data_(Bitwise(data, nbits/8, false)), nhf_(nhf), nkeys_(nkeys), seeds_(seeds), nmod_(nbits) {}

        void init();
        bool AddKeys(const vector<Bitwise> &keys, const vector<uint16_t> & runids = vector<uint16_t>());
        bool AddKeys_len(const vector<Bitwise> &keys, const vector<uint16_t> & runids = vector<uint16_t>());
        bool Query(const Bitwise &key);
        bool Query_len(const Bitwise &key);
        void printStats() const {
            assert(keep_stats);
            printf("#queries: %lu, #positives: %lu\n", nqueries_, npositives_);
        }
        size_t mem() const{
            return sizeof(*this) + sizeof(size_t) * seeds_.size() + data_.size()/8;
        }
};

#ifdef USE_DTL
class DtlBlockedBloomFilter final: public Filter {
    private:
        Bitwise data_;
        size_t nkeys_;
//        dtl::blocked_bloomfilter<uint64_t> *b_;
        dtl::bbf_64 *b_;

    public:
        DtlBlockedBloomFilter(size_t nbits): data_(Bitwise(false, ((nbits+7)/8)*8)) {}
        //~DtlBlockedBloomFilter() { delete b_; }

        bool AddKeys(const vector<Bitwise> &keys, const vector<uint16_t> & runids);
        bool Query(const Bitwise &key);
};
#endif


template<class FilterClass, bool keep_stats=false>
class Rosetta final: public Filter {
    private:
        size_t diffidence_, maxlen_, nkeys_, diffidence_level_;
        vector<FilterClass*> bfs_;
        function<vector<size_t> (vector<size_t>)> get_nbits_;

    public:
        size_t nqueries_ = 0;
        size_t npositives_ = 0;
        vector<size_t> qdist_;

        Rosetta(size_t diffidence, size_t diffidence_level, function<vector<size_t> (vector<size_t>)> get_nbits): diffidence_(diffidence), maxlen_(0), nkeys_(0), diffidence_level_(diffidence_level), get_nbits_(get_nbits) {
            static_assert(is_base_of<Filter, FilterClass>::value, "DST template argument must be a filter");
        }
        ~Rosetta(){
            for(auto &bf: bfs_)
                delete bf;
        }

        bool AddKeys(const vector<Bitwise> &keys, const vector<uint16_t> & runids);
        bool Doubt(Bitwise *idx, size_t &C, size_t level, size_t maxlevel);
        Bitwise *GetFirst(const Bitwise &from, const Bitwise &to);
        Bitwise *Seek(const Bitwise &from);
        bool Query(const Bitwise &key);
        bool Query(const Bitwise &from, const Bitwise &to);
        void printStats() const {
            assert(keep_stats);
            printf("DST total stats: #queries: %lu, #positives: %lu\n", nqueries_, npositives_);
            printf("Stats for each bf:\n");
            for (auto &bf: bfs_) {
                printf("\t");
                bf->printStats();
            }
            printf("Query distribution:\n");
            for (auto &i: qdist_) {
                printf("\t%lu\n", i);
            }
        }
        size_t mem() const;
};

namespace vacuum{

#define memcle(a) memset(a, 0, sizeof(a))
#define sqr(a) ((a) * (a))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ROUNDDOWN(a, b) ((a) - ((a) % (b)))
#define ROUNDUP(a, b) ROUNDDOWN((a) + (b - 1), b)


template <typename fp_t, int fp_len>
class vFilter : public Filter {
public:
    long long n;  // number of buckets
    int m;        // number of slots per bucket
    uint64_t memory_consumption;
    virtual void init(int _n, int _m, int _max_kick_steps) = 0;
    virtual void clear() = 0;
    virtual bool insert(uint64_t ele) = 0;
    virtual bool lookup(uint64_t ele) = 0;
    virtual bool del(uint64_t ele) = 0;
    uint64_t position_hash(uint64_t ele);  // hash to range [0, n - 1]
    virtual double get_load_factor() { return 0; }
    virtual double get_full_bucket_factor() { return 0; }
    virtual void debug_test() {}
};

template <typename fp_t, int fp_len>
class SemiSortCuckooFilter : public vFilter<fp_t, fp_len> {
public:
    int max_2_power;
    int big_seg;
    int len[4];
    virtual void init(int _n, int _m, int _max_kick_steps);
    void clear();
    virtual bool insert(uint64_t ele);
    bool insert(uint64_t ele, uint16_t runid);
    bool lookup(uint64_t ele);
    bool lookupIO(uint64_t ele, const Bitwise &key);
    virtual bool del(uint64_t ele);
    double get_load_factor();
    double get_full_bucket_factor();
    double get_bits_per_item();

    bool debug_flag = false;
    bool balance = true;
    uint32_t* T;
    static uint32_t encode_table[1 << 16];
    static uint32_t decode_table[1 << 16];

    uint16_t* rid; // for simulation

    ~SemiSortCuckooFilter() { delete T; delete rid; }

    int filled_cell;
    int full_bucket;
    int max_kick_steps;

    function<bool (uint16_t, const Bitwise&)>* io_sim_;

    fp_t fingerprint(uint64_t ele);  // 32-bit to 'fp_len'-bit fingerprint

    // interface for semi-sorted bucket
    void get_bucket(int pos, fp_t* store);
    void set_bucket(int pos, fp_t* sotre);
    void test_bucket();
    void make_balance();
    inline int high_bit(fp_t fp);
    inline int low_bit(fp_t fp);
    inline void sort_pair(fp_t& a, fp_t& b);

    virtual int alternate(int pos, fp_t fp) = 0;  // get alternate position
    virtual int insert_to_bucket(
        fp_t* store, fp_t fp);  // insert one fingerprint to bucket [pos]
    virtual int lookup_in_bucket(
        int pos, fp_t fp);  // lookup one fingerprint in  bucket [pos]
    virtual int lookupIO_in_bucket(
        int pos, fp_t fp, const Bitwise& key);
    virtual int del_in_bucket(
        int pos, fp_t fp);  // lookup one fingerprint in  bucket [pos]
};

int upperpower2(int x);

// solve equation : 1 + x(logc - logx + 1) - c = 0
double F_d(double x, double c);
double F(double x, double c);
double solve_equation(double c);
double balls_in_bins_max_load(double balls, double bins);

int proper_alt_range(int M, int i, int* len);

template <typename fp_t, int fp_len>
class VacuumFilter final : public SemiSortCuckooFilter<fp_t, fp_len> {
private:
    virtual int alternate(int pos, fp_t fp);
    size_t seed;
public:
    VacuumFilter(size_t nbits);
    VacuumFilter(size_t nbits, function<bool (uint16_t, const Bitwise&)>* io_sim);
    bool AddKeys(const vector<Bitwise> &keys, const vector<uint16_t> & runids);
    bool AddKeys_len(const vector<Bitwise> &keys, const vector<uint16_t> & runids);
    bool Query(const Bitwise &key);
    bool QueryIO(const Bitwise &key);
    bool Query_len(const Bitwise &key);
    bool QueryIO_len(const Bitwise &key);
    size_t mem() const;

    void printStats() const {}
};

}

template<class FilterClass>
class SplittedRosetta final: public Filter {
    private:
        size_t maxlen_, nkeys_;
        vector<BloomFilter<>*> bfs_;
        vector<FilterClass*> cks_;
        function<pair<vector<size_t>, vector<size_t>> (vector<size_t>, vector<size_t>, size_t, size_t, uint64_t)> get_nbits_;
        function<bool (uint16_t, const Bitwise&)> io_sim_;
        uint64_t ck_mask_;
        size_t bf_max_, ck_max_;

    public:

        SplittedRosetta(size_t maxlen, size_t bf_max, size_t ck_max, uint64_t ck_mask, function<pair<vector<size_t>, vector<size_t>> (vector<size_t>, vector<size_t>, size_t, size_t, uint64_t)> get_nbits, function<bool (uint16_t, const Bitwise&)> io_sim):
        maxlen_(maxlen), nkeys_(0), bf_max_(bf_max), ck_max_(ck_max), ck_mask_(ck_mask), get_nbits_(get_nbits), io_sim_(io_sim) {
            static_assert(is_base_of<Filter, FilterClass>::value, "DST template argument must be a filter");
        }
        ~SplittedRosetta(){
            for(auto &bf: bfs_)
                delete bf;
            for(auto &ck: cks_)
                delete ck;
        }

        bool AddKeys(const vector<Bitwise> &keys, const vector<uint16_t> & runids);
        bool Doubt(Bitwise *idx, size_t level);
        Bitwise *GetFirst(const Bitwise &from, const Bitwise &to) {}
        Bitwise *Seek(const Bitwise &from);
        bool Query(const Bitwise &key) {}
        bool Query(const Bitwise &from, const Bitwise &to) {}
        size_t mem() const;
};

#endif
