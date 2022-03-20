#include <vector>
#include <cassert>
#include <memory>
#include <cstring>
#include <cmath>
#include <random>

#include "dst.h"

vector<size_t> calc_dst(vector<size_t> dist, double bpk, vector<size_t> qdist, size_t cutoff) {
    vector<size_t> out(dist.size(), 0);
    size_t bloom_cost = (sizeof(BloomFilter<>) + sizeof(size_t))*8;
    size_t totbits = (size_t)((long double)bpk*dist[dist.size()-1]);

    long double c=0;
    size_t mqdist = 1000000000;
    for (size_t i=0; i<qdist.size(); ++i)
        if (qdist[i] > 0) {
            totbits = totbits>bloom_cost?totbits-bloom_cost:0;
            mqdist = min(mqdist, qdist[i]);
        }

    long double lo=1/(long double)mqdist, hi=100000, mid;
    size_t bitsused = 0;
    while (fabs(hi-lo) > 1/(long double)mqdist/10) {
        mid = (lo+hi)/2;

        bitsused = 0;
        bool ok = true;
        for (size_t i=0; i<dist.size()-1-cutoff; ++i)
            if (dist[i] < (1ULL<<(i+1)) and qdist[i] > 0) {
                long double next_fpr = 1/(mid*qdist[i+1]);
                long double this_fpr = 1/(mid*qdist[i]);
                long double bloom_fpr = this_fpr/(2-next_fpr)/next_fpr;
                bloom_fpr = min(bloom_fpr, (long double)1);
                bloom_fpr = max(bloom_fpr, (long double)1e-12);
                bitsused += (size_t)ceil((long double)dist[i]*1.44*log2(1/bloom_fpr));
                if (bitsused > totbits) {
                    ok = false;
                    break;
                }
            }
        bitsused += (size_t)ceil((long double)dist[dist.size()-1]*1.44*log2(mid*qdist[dist.size()-1]));
        if (ok and bitsused <= totbits)
            lo = mid;
        else
            hi = mid;
    }
    if (lo < 1e-12) {
        printf("Not enough memory allowance\n");
        exit(0);
    }
    c = lo;

    bitsused = 0;
    for (size_t i=0; i<dist.size()-1-cutoff; ++i)
        if (dist[i] < (1ULL<<(i+1)) and qdist[i] > 0 and bitsused < totbits) {
            long double next_fpr = 1/(c*qdist[i+1]);
            long double this_fpr = 1/(c*qdist[i]);
            long double bloom_fpr = this_fpr/(2-next_fpr)/next_fpr;
            bloom_fpr = min(bloom_fpr, (long double)1);
            bloom_fpr = max(bloom_fpr, (long double)1e-12);
            out[i] = min((size_t)ceil((long double)dist[i]*1.44*log2(1/bloom_fpr)), totbits-bitsused);
            bitsused += out[i];
        }
    out[dist.size()-1] = bitsused<=totbits?totbits-bitsused:0;
    bitsused += out[dist.size()-1];
    if (bitsused > totbits) {
        printf("Not enough memory allowance, bitsused: %lu, totbits: %lu\n", bitsused, totbits);
        exit(0);
    }
    return out;
}

template<bool keep_stats>
void BloomFilter<keep_stats>::AddKeys(const vector<Bitwise> &keys) {
    if (data_.size() == 0) return;
    nkeys_ = keys.size();
    nhf_ = (size_t)(round(log(2)*data_.size()/nkeys_));
    nhf_ = (nhf_==0?1:nhf_);
    seeds_.resize(nhf_);
    mt19937 gen(1337);
    for (size_t i=0; i<nhf_; ++i)
        seeds_[i] = gen();

    for (auto &key: keys)
        for (size_t i=0; i<nhf_; ++i)
            data_.set(key.hash(seeds_[i], nmod_), 1);
}

template<bool keep_stats>
bool BloomFilter<keep_stats>::Query(const Bitwise &key) {
    if (data_.size() == 0) return true;
    bool out=true;
    for (size_t i=0; i<nhf_ && out; ++i)
        out &= data_.get(key.hash(seeds_[i], nmod_));
    if (keep_stats) {
        nqueries_ += 1;
        npositives_ += out;
    }
    return out;
}

template<bool keep_stats>
pair<uint8_t*, size_t> BloomFilter<keep_stats>::serialize() const {
    size_t len = sizeof(size_t) + sizeof(size_t) + sizeof(size_t) + seeds_.size()*sizeof(size_t) + (data_.size()+7)/8;
    uint8_t *out = new uint8_t[len];
    uint8_t *pos = out;

    memcpy(pos, &nkeys_, sizeof(size_t));
    pos += sizeof(size_t);
    memcpy(pos, &nhf_, sizeof(size_t));
    pos += sizeof(size_t);
    size_t tmp = data_.size();
    memcpy(pos, &tmp, sizeof(size_t));
    pos += sizeof(size_t);
    memcpy(pos, seeds_.data(), seeds_.size()*sizeof(size_t));
    pos += seeds_.size()*sizeof(size_t);
    memcpy(pos, data_.data(), (data_.size()+7)/8);
    return {out, len};
}

template<bool keep_stats>
pair<BloomFilter<keep_stats>*, size_t> BloomFilter<keep_stats>::deserialize(uint8_t* ser) {
    uint8_t *pos = ser;

    size_t nkeys = ((size_t*)pos)[0];
    size_t nhf = ((size_t*)pos)[1];
    size_t nbits = ((size_t*)pos)[2];
    pos += 3*sizeof(size_t);
    vector<size_t> seeds((size_t*)pos, (size_t*)pos+nhf);
    pos += nhf*sizeof(size_t);

    size_t len = pos-ser + nbits/8;
    return {new BloomFilter(nbits, pos, nhf, nkeys, seeds), len};
}

#ifdef USE_DTL
void DtlBlockedBloomFilter::AddKeys(const vector<Bitwise> &keys) {
    if (data_.size() == 0) return;
    nkeys_ = keys.size();
    size_t k = data_.size()/nkeys_;
    if (k == 0) {
        k = 1;
    }
    dtl::bbf_64::force_unroll_factor(0);
    b_ = new dtl::bbf_64(data_.size(), k, 1, 1);
//    b_ = new dtl::blocked_bloomfilter<uint64_t>(data_.size(), k, 4, 1, dtl::blocked_bloomfilter_tune());
    for (auto &key: keys) {
        b_->insert((uint64_t*)data_.data(), key.to_uint64());
    }
}

bool DtlBlockedBloomFilter::Query(const Bitwise &key) {
    if (data_.size() == 0) return true;
    return b_->contains((uint64_t*)data_.data(), key.to_uint64());
}

pair<uint8_t*, size_t> DtlBlockedBloomFilter::serialize() const {
    size_t len = sizeof(size_t) + sizeof(size_t) + (data_.size()+7)/8;
    uint8_t *out = new uint8_t[len];
    uint8_t *pos = out;

    memcpy(pos, &nkeys_, sizeof(size_t));
    pos += sizeof(size_t);
    size_t tmp = data_.size();
    memcpy(pos, &tmp, sizeof(size_t));
    pos += sizeof(size_t);
    memcpy(pos, data_.data(), (data_.size()+7)/8);
    return {out, len};
}
#endif

template<class FilterClass, bool keep_stats>
void Rosetta<FilterClass, keep_stats>::AddKeys(const vector<Bitwise> &keys) {
    nkeys_ = keys.size();
    vector<size_t> distribution;
    vector<vector<Bitwise>> bloom_keys;
    maxlen_ = 0;
    for (size_t i=0; i<nkeys_; ++i)
        maxlen_ = max(maxlen_, keys[i].size());
    distribution.resize(maxlen_, 0);
    bloom_keys.resize(maxlen_);
    for (size_t i=0; i<nkeys_; ++i) {
        size_t lcp = i>0?keys[i-1].lcp(keys[i]):0;
        for (size_t j=lcp; j<maxlen_; ++j)
            ++distribution[j];
    }
    vector<size_t> nbits = get_nbits_(distribution);

    for (size_t j=0; j<maxlen_; ++j)
        if (nbits[j] > 0)
            bloom_keys[j].reserve(distribution[j]);
    for (size_t i=0; i<nkeys_; ++i) {
        size_t lcp = i>0?keys[i-1].lcp(keys[i]):0;
        for (size_t j=lcp; j<maxlen_; ++j)
            if (nbits[j] > 0)
                bloom_keys[j].emplace_back(Bitwise(keys[i], j+1));
    }

    if (keep_stats)
        qdist_.resize(maxlen_, 0);

    bfs_.reserve(maxlen_);
    for (size_t i=0; i<maxlen_; ++i) {
        bfs_.emplace_back(new FilterClass(nbits[i]));
        bfs_[i]->AddKeys(bloom_keys[i]);
    }
}

template<class FilterClass, bool keep_stats>
bool Rosetta<FilterClass, keep_stats>::Doubt(Bitwise *idx, size_t &C, size_t level, size_t maxlevel) {
    if (level >= maxlen_)
        return false;
    if (not bfs_[level]->Query(Bitwise(*idx, level+1)))
        return false;
    --C;
    if (level == maxlevel-1 or C == 0) {
        C = 0;
        return true;
    }
    bool old = idx->get(level+1);
    idx->set(level+1, 0);
    if (Doubt(idx, C, level+1, maxlevel))
        return true;
    idx->set(level+1, 1);
    if (Doubt(idx, C, level+1, maxlevel))
        return true;
    idx->set(level+1, old);
    return false;
}

template<class FilterClass, bool keep_stats>
Bitwise *Rosetta<FilterClass, keep_stats>::GetFirst(const Bitwise &from, const Bitwise &to) {
    Bitwise tfrom(from, min(from.size(), maxlen_)),
            tto(to, min(to.size(), maxlen_));
    size_t lcp = tfrom.lcp(tto);
    size_t C = diffidence_;

    if (!keep_stats) {
        Bitwise *tmp = new Bitwise(tto.data(), tto.size()/8);
        for (size_t i=0; i<lcp; ++i) {
            if (!Doubt(tmp, C, i, i+1))
                return tmp;
            memcpy(tmp->data(), tto.data(), tto.size()/8);
        }
        delete tmp;
    }

    bool carry = false;
    Bitwise *out;
    if (tfrom.size() > tto.size())
        out = new Bitwise(tfrom.data(), tfrom.size()/8);
    else {
        out = new Bitwise(false, tto.size());
        memcpy(out->data(), tfrom.data(), tfrom.size()/8);
    }
    if (lcp == maxlen_)
        return out;
    for (size_t i=tfrom.size()-1; i>lcp; --i) {
        if (carry) {
            if (out->get(i) == 0)
                carry = false;
            out->flip(i);
        }
        if (out->get(i) == 1) {
            if (keep_stats)
                ++qdist_[i];
            if (Doubt(out, C, i, min(min(maxlen_, out->size()), i+diffidence_level_)))
                return out;
            out->set(i, 0);
            carry = true;
        }
    }

    if (keep_stats and !carry and lcp < maxlen_)
        ++qdist_[lcp];
    if (!carry and Doubt(out, C, lcp, min(min(maxlen_, out->size()), lcp+diffidence_level_)))
        return out;
    out->set(lcp, 1);

    for (size_t i=lcp+1; i<tto.size(); ++i)
        if (tto.get(i) == 1) {
            if (keep_stats)
                ++qdist_[i];
            if (Doubt(out, C, i, min(min(maxlen_, out->size()), i+diffidence_level_)))
                return out;
            out->set(i, 1);
        }
    return out;
}

template<class FilterClass, bool keep_stats>
Bitwise *Rosetta<FilterClass, keep_stats>::Seek(const Bitwise &from) {
    Bitwise tfrom(from, min(from.size(), maxlen_));
    size_t C = diffidence_;
    bool carry = false;
    Bitwise *out = new Bitwise(tfrom.data(), tfrom.size()/8);
    for (size_t i=tfrom.size()-1; i>=0; --i) {
        if (carry) {
            if (out->get(i) == 0)
                carry = false;
            out->flip(i);
        }
        if (out->get(i) == 1) {
            if (keep_stats)
                ++qdist_[i];
            if (Doubt(out, C, i, min(min(maxlen_, out->size()), i+diffidence_level_)))
                return out;
            out->set(i, 0);
            carry = true;
        }
    }
    if (!carry){
        if (keep_stats)
            qdist_[0] += 2;
        if(Doubt(out, C, 0, min(min(maxlen_, out->size()), diffidence_level_)))
            return out;
        out->set(0, 1);
        if(Doubt(out, C, 0, min(min(maxlen_, out->size()), diffidence_level_)))
            return out;
    }
    return out;
}

template<class FilterClass, bool keep_stats>
bool Rosetta<FilterClass, keep_stats>::Query(const Bitwise &from, const Bitwise &to) {
    Bitwise *qry = GetFirst(from, to);
    size_t lcp = to.lcp(*qry);
    delete qry;
    if (keep_stats)
        nqueries_ += 1;
    if (lcp < min(to.size(), maxlen_)) {
        if (keep_stats)
            npositives_ += 1;
        return true;
    }
    return false;
}

template<class FilterClass, bool keep_stats>
bool Rosetta<FilterClass, keep_stats>::Query(const Bitwise &key){
    if (keep_stats)
        ++qdist_[maxlen_-1];
    return bfs_[maxlen_-1]->Query(Bitwise(key, maxlen_));
}

template<class FilterClass, bool keep_stats>
pair<uint8_t*, size_t> Rosetta<FilterClass, keep_stats>::serialize() const {
    vector<pair<uint8_t*, size_t>> serial_Bloom;
    size_t len = 4*sizeof(size_t);
    for (auto &bf: bfs_) {
        auto ser = bf->serialize();
        len += ser.second;
        serial_Bloom.emplace_back(ser);
    }
    uint8_t *out = new uint8_t[len];
    uint8_t *pos = out;
    memcpy(pos, &diffidence_, sizeof(size_t));
    pos += sizeof(size_t);
    memcpy(pos, &diffidence_level_, sizeof(size_t));
    pos += sizeof(size_t);
    memcpy(pos, &maxlen_, sizeof(size_t));
    pos += sizeof(size_t);
    memcpy(pos, &nkeys_, sizeof(size_t));
    pos += sizeof(size_t);
    for (auto &ser: serial_Bloom) {
        memcpy(pos, ser.first, ser.second);
        delete[] ser.first;
        pos += ser.second;
    }
    return {out, len};
}
template<class FilterClass, bool keep_stats>
pair<Rosetta<FilterClass, keep_stats>*, size_t> Rosetta<FilterClass, keep_stats>::deserialize(uint8_t* ser) {
    uint8_t* pos = ser;
    size_t diffidence = ((size_t*)pos)[0];
    size_t diffidence_level = ((size_t*)pos)[1];
    Rosetta<FilterClass, keep_stats>* out = new Rosetta<FilterClass, keep_stats>(diffidence, diffidence_level, [](vector<size_t> distribution) -> vector<size_t> { assert(false); return distribution;});
    out->maxlen_ = ((size_t*)pos)[2];
    if (keep_stats)
        out->qdist_.resize(out->maxlen_, 0);
    out->nkeys_ = ((size_t*)pos)[3];
    pos += 4*sizeof(size_t);
    for (size_t i=0; i<out->maxlen_; ++i) {
        pair<FilterClass*,size_t> tmp = FilterClass::deserialize(pos);
        out->bfs_.emplace_back(tmp.first);
        delete tmp.first;
        pos += tmp.second;
    }
    return {out, pos-ser};
}

template<class FilterClass, bool keep_stats>
size_t Rosetta<FilterClass, keep_stats>::mem() const{
    size_t s = sizeof(*this);
    for(auto &b: bfs_)
        s += b->mem();
    return s;
}


template<class FilterClass, bool keep_stats>
void UnlayeredRosetta<FilterClass, keep_stats>::AddKeys(const vector<Bitwise> &keys) {
    nkeys_ = keys.size();
    vector<size_t> distribution;
    vector<vector<Bitwise>> bloom_keys;
    maxlen_ = 0;
    for (size_t i=0; i<nkeys_; ++i)
        maxlen_ = max(maxlen_, keys[i].size());
    distribution.resize(maxlen_, 0);
    bloom_keys.resize(maxlen_);
    for (size_t i=0; i<nkeys_; ++i) {
        size_t lcp = i>0?keys[i-1].lcp(keys[i]):0;
        for (size_t j=lcp; j<maxlen_; ++j)
            ++distribution[j];
    }
    size_t nbits = get_nbits_(distribution);

    for (size_t j=0; j<maxlen_; ++j)
        bloom_keys[j].reserve(distribution[j]);
    for (size_t i=0; i<nkeys_; ++i) {
        size_t lcp = i>0?keys[i-1].lcp(keys[i]):0;
        for (size_t j=lcp; j<maxlen_; ++j)
            bloom_keys[j].emplace_back(Bitwise(keys[i], j+1));
    }

    if (keep_stats)
        qdist_.resize(maxlen_, 0);

    bfs_ = new FilterClass(nbits);
    for (size_t i=0; i<maxlen_; ++i)
        bfs_->AddKeys_len(bloom_keys[i]);
}

template<class FilterClass, bool keep_stats>
bool UnlayeredRosetta<FilterClass, keep_stats>::Doubt(Bitwise *idx, size_t &C, size_t level, size_t maxlevel) {
    if (level >= maxlen_)
        return false;
    if (not bfs_->Query_len(Bitwise(*idx, level+1)))
        return false;
    --C;
    if (level == maxlevel-1 or C == 0) {
        C = 0;
        return true;
    }
    bool old = idx->get(level+1);
    idx->set(level+1, 0);
    if (Doubt(idx, C, level+1, maxlevel))
        return true;
    idx->set(level+1, 1);
    if (Doubt(idx, C, level+1, maxlevel))
        return true;
    idx->set(level+1, old);
    return false;
}

template<class FilterClass, bool keep_stats>
Bitwise *UnlayeredRosetta<FilterClass, keep_stats>::GetFirst(const Bitwise &from, const Bitwise &to) {
    Bitwise tfrom(from, min(from.size(), maxlen_)),
            tto(to, min(to.size(), maxlen_));
    size_t lcp = tfrom.lcp(tto);
    size_t C = diffidence_;

    if (!keep_stats) {
        Bitwise *tmp = new Bitwise(tto.data(), tto.size()/8);
        for (size_t i=0; i<lcp; ++i) {
            if (!Doubt(tmp, C, i, i+1))
                return tmp;
            memcpy(tmp->data(), tto.data(), tto.size()/8);
        }
        delete tmp;
    }

    bool carry = false;
    Bitwise *out;
    if (tfrom.size() > tto.size())
        out = new Bitwise(tfrom.data(), tfrom.size()/8);
    else {
        out = new Bitwise(false, tto.size());
        memcpy(out->data(), tfrom.data(), tfrom.size()/8);
    }
    if (lcp == maxlen_)
        return out;
    for (size_t i=tfrom.size()-1; i>lcp; --i) {
        if (carry) {
            if (out->get(i) == 0)
                carry = false;
            out->flip(i);
        }
        if (out->get(i) == 1) {
            if (keep_stats)
                ++qdist_[i];
            if (Doubt(out, C, i, min(min(maxlen_, out->size()), i+diffidence_level_)))
                return out;
            out->set(i, 0);
            carry = true;
        }
    }

    if (keep_stats and !carry and lcp < maxlen_)
        ++qdist_[lcp];
    if (!carry and Doubt(out, C, lcp, min(min(maxlen_, out->size()), lcp+diffidence_level_)))
        return out;
    out->set(lcp, 1);

    for (size_t i=lcp+1; i<tto.size(); ++i)
        if (tto.get(i) == 1) {
            if (keep_stats)
                ++qdist_[i];
            if (Doubt(out, C, i, min(min(maxlen_, out->size()), i+diffidence_level_)))
                return out;
            out->set(i, 1);
        }
    return out;
}

template<class FilterClass, bool keep_stats>
Bitwise *UnlayeredRosetta<FilterClass, keep_stats>::Seek(const Bitwise &from) {
    Bitwise tfrom(from, min(from.size(), maxlen_));
    size_t C = diffidence_;
    bool carry = false;
    Bitwise *out = new Bitwise(tfrom.data(), tfrom.size()/8);
    for (size_t i=tfrom.size()-1; i>=0; --i) {
        if (carry) {
            if (out->get(i) == 0)
                carry = false;
            out->flip(i);
        }
        if (out->get(i) == 1) {
            if (keep_stats)
                ++qdist_[i];
            if (Doubt(out, C, i, min(min(maxlen_, out->size()), i+diffidence_level_)))
                return out;
            out->set(i, 0);
            carry = true;
        }
    }
    if (!carry){
        if (keep_stats)
            qdist_[0] += 2;
        if(Doubt(out, C, 0, min(min(maxlen_, out->size()), diffidence_level_)))
            return out;
        out->set(0, 1);
        if(Doubt(out, C, 0, min(min(maxlen_, out->size()), diffidence_level_)))
            return out;
    }
    return out;
}

template<class FilterClass, bool keep_stats>
bool UnlayeredRosetta<FilterClass, keep_stats>::Query(const Bitwise &from, const Bitwise &to) {
    Bitwise *qry = GetFirst(from, to);
    size_t lcp = to.lcp(*qry);
    delete qry;
    if (keep_stats)
        nqueries_ += 1;
    if (lcp < min(to.size(), maxlen_)) {
        if (keep_stats)
            npositives_ += 1;
        return true;
    }
    return false;
}

template<class FilterClass, bool keep_stats>
bool UnlayeredRosetta<FilterClass, keep_stats>::Query(const Bitwise &key){
    if (keep_stats)
        ++qdist_[maxlen_-1];
    return bfs_->Query_len(Bitwise(key, maxlen_));
}

template<class FilterClass, bool keep_stats>
size_t UnlayeredRosetta<FilterClass, keep_stats>::mem() const{
    return sizeof(*this) + bfs_->mem();
}


namespace vacuum{

template <typename fp_t, int fp_len>
uint64_t vFilter<fp_t, fp_len>::position_hash(uint64_t ele) {
    return (ele>>32) % n;
}

int upperpower2(int x) {
    int ret = 1;
    for (; ret < x;) ret <<= 1;
    return ret;
}

// solve equation : 1 + x(logc - logx + 1) - c = 0
double F_d(double x, double c) { return log(c) - log(x); }
double F(double x, double c) { return 1 + x * (log(c) - log(x) + 1) - c; }
double solve_equation(double c) {
    double x = c + 0.1;
    while (abs(F(x, c)) > 0.001) x -= F(x, c) / F_d(x, c);
    return x;
}
double balls_in_bins_max_load(double balls, double bins) {
    double m = balls;
    double n = bins;
    if (n == 1) return m;

    double c = m / (n * log(n));
    // A more accurate bound..
    if (c < 5) {
        double dc = solve_equation(c);
        double ret = (dc - 1 + 2) * log(n);
        return ret;
    }

    double ret = (m / n) + 1.5 * sqrt(2 * m / n * log(n));
    return ret;
}

int proper_alt_range(int M, int i, int* len) {
    double b = 4;      // slots per bucket
    double lf = 0.95;  // target load factor
    int alt_range = 8;
    for (; alt_range < M;) {
        double f = (4 - i) * 0.25;
        if (balls_in_bins_max_load(f * b * lf * M, M * 1.0 / alt_range) <
            0.97 * b * alt_range) {
            return alt_range;
        }
        alt_range <<= 1;
    }
    return alt_range;
}
template<> uint32_t SemiSortCuckooFilter<uint16_t, 16>::encode_table[1 << 16] = {};
template<> uint32_t SemiSortCuckooFilter<uint16_t, 16>::decode_table[1 << 16] = {};
template<> uint32_t SemiSortCuckooFilter<uint8_t, 8>::encode_table[1 << 16] = {};
template<> uint32_t SemiSortCuckooFilter<uint8_t, 8>::decode_table[1 << 16] = {};
//+++ We need to make the 'max_item' to be the number of buckets
template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::init(int max_item, int _m, int _step) {
//    int _n = (max_item / 0.96 / 4);
    int _n = MAX((max_item / 4), 1);

    if (false && _n < 10000) {
        if (_n < 256)
            big_seg = (upperpower2(_n));
        else
            big_seg = (upperpower2(_n / 4));
        _n = ROUNDUP(_n, big_seg);
        len[0] = big_seg - 1;
        len[1] = big_seg - 1;
        len[2] = big_seg - 1;
        len[3] = big_seg - 1;
    } else {
        big_seg = 0;
        big_seg = max(big_seg, proper_alt_range(_n, 0, len));
        int new_n = ROUNDUP(_n, big_seg);
        _n = new_n;

        big_seg--;
        len[0] = big_seg;
        for (int i = 1; i < 4; i++) len[i] = proper_alt_range(_n, i, len) - 1;
        len[3] = (len[3] + 1) * 2 - 1;
    }

    this->n = _n;
    this->m = _m;
    this->max_kick_steps = _step;
    this->filled_cell = 0;
    this->full_bucket = 0;

    uint64_t how_many_bit = (uint64_t)this->n * this->m * (fp_len - 1);
    this->memory_consumption =
        ROUNDUP(how_many_bit + 64, 8) / 8 + 8;  // how many bytes !

    max_2_power = 1;
    for (; max_2_power * 2 < _n;) max_2_power <<= 1;
    this->T = new uint32_t[this->memory_consumption]();

    int index = 0;
    for (int i = 0; i < 16; i++)
        for (int j = 0; j < ((i == 0) ? 1 : i + 1); j++)
            for (int k = 0; k < ((j == 0) ? 1 : j + 1); k++)
                for (int l = 0; l < ((k == 0) ? 1 : k + 1); l++) {
                    int plain_bit = (i << 12) + (j << 8) + (k << 4) + l;
                    encode_table[plain_bit] = index;
                    decode_table[index] = plain_bit;
                    ++index;
                }
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::clear() {
    this->filled_cell = 0;
    memset(this->T, 0, this->memory_consumption);
}

template <typename fp_t, int fp_len>
fp_t SemiSortCuckooFilter<fp_t, fp_len>::fingerprint(uint64_t ele) {
    fp_t h = (ele & 0x0000ffffu) % ((1ull << fp_len) - 1) + 1;
    return h;
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::get_bucket(int pos, fp_t* store) {
    // Default :
    //
    // Little Endian Store
    // Store by uint32_t
    // store[this -> m] = bucket number

    // 1. read the endcoded bits from memory

    int bucket_length = (fp_len - 1) * 4;
    uint64_t start_bit_pos = (uint64_t)pos * bucket_length;
    uint64_t end_bit_pos = start_bit_pos + bucket_length - 1;
    uint64_t result = 0;

    if (ROUNDDOWN(start_bit_pos, 64) == ROUNDDOWN(end_bit_pos, 64)) {
        uint64_t unit = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64];
        int reading_lower_bound = start_bit_pos & 63;
        int reading_upper_bound = end_bit_pos & 63;

        result = ((uint64_t)unit & ((-1ULL) >> (63 - reading_upper_bound))) >>
                 reading_lower_bound;
    } else {
        uint64_t unit1 = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64];
        uint64_t unit2 = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64 + 1];

        int reading_lower_bound = start_bit_pos & 63;
        int reading_upper_bound = end_bit_pos & 63;

        uint64_t t1 = unit1 >> reading_lower_bound;
        uint64_t t2 = (unit2 & ((1ULL << (reading_upper_bound + 1)) - 1))
                      << (64 - reading_lower_bound);
        result = t1 + t2;
    }

    // 2. read the 4 elements from the encoded bits
    // We use 12 bits to store the 16 most significant bits for the items in
    // bucket, 4 bits per item the low bits are stored in the remaining bits
    //
    // For example, 8 bits per item , require 28 bits to store:
    //
    // Original :
    //
    // hhhh llll
    // hhhh llll
    // hhhh llll
    // hhhh llll
    //
    // encoded :
    //
    //
    // 0 - 11                       12 - 15    16 - 19  20-23   24 - 27
    // HHHHHHHHHHHH                 llll       llll     llll    llll
    //  encoded high bit(12 bits)   item 0     item 1   item 2  item 3
    //
    int decode_result = decode_table[result >> (4 * (fp_len - 4))];

    store[3] = (result & ((1 << (fp_len - 4)) - 1)) +
               ((decode_result & ((1 << 4) - 1)) << (fp_len - 4));
    store[2] = ((result >> (1 * (fp_len - 4))) & ((1 << (fp_len - 4)) - 1)) +
               (((decode_result >> 4) & ((1 << 4) - 1)) << (fp_len - 4));
    store[1] = ((result >> (2 * (fp_len - 4))) & ((1 << (fp_len - 4)) - 1)) +
               (((decode_result >> 8) & ((1 << 4) - 1)) << (fp_len - 4));
    store[0] = ((result >> (3 * (fp_len - 4))) & ((1 << (fp_len - 4)) - 1)) +
               (((decode_result >> 12) & ((1 << 4) - 1)) << (fp_len - 4));

    store[4] = 0;
    store[4] += store[0] != 0;
    store[4] += store[1] != 0;
    store[4] += store[2] != 0;
    store[4] += store[3] != 0;
}

template <typename fp_t, int fp_len>
inline void SemiSortCuckooFilter<fp_t, fp_len>::sort_pair(fp_t& a, fp_t& b) {
    if ((a) < (b)) swap(a, b);
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::set_bucket(int pos, fp_t* store) {
    // 0. sort store ! descendant order >>>>>>

    sort_pair(store[0], store[2]);
    sort_pair(store[1], store[3]);
    sort_pair(store[0], store[1]);
    sort_pair(store[2], store[3]);
    sort_pair(store[1], store[2]);

    // 1. compute the encode

    uint64_t high_bit = 0;
    uint64_t low_bit = 0;

    low_bit =
        (store[3] & ((1 << (fp_len - 4)) - 1)) +
        ((store[2] & ((1 << (fp_len - 4)) - 1)) << (1 * (fp_len - 4))) +
        (((uint64_t)store[1] & ((1 << (fp_len - 4)) - 1)) << (2 * (fp_len - 4))) +
        (((uint64_t)store[0] & ((1 << (fp_len - 4)) - 1)) << (3 * (fp_len - 4)));

    high_bit = ((store[3] >> (fp_len - 4)) & ((1 << 4) - 1)) +
               (((store[2] >> (fp_len - 4)) & ((1 << 4) - 1)) << 4) +
               (((store[1] >> (fp_len - 4)) & ((1 << 4) - 1)) << 8) +
               (((store[0] >> (fp_len - 4)) & ((1 << 4) - 1)) << 12);


    // 2. store into memory
    uint64_t high_encode = encode_table[high_bit];
    uint64_t all_encode = (high_encode << (4 * (fp_len - 4))) + low_bit;

    int bucket_length = (fp_len - 1) * 4;
    uint64_t start_bit_pos = (uint64_t)pos * bucket_length;
    uint64_t end_bit_pos = start_bit_pos + bucket_length - 1;

    if (ROUNDDOWN(start_bit_pos, 64) == ROUNDDOWN(end_bit_pos, 64)) {
        uint64_t unit = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64];
        int writing_lower_bound = start_bit_pos & 63;
        int writing_upper_bound = end_bit_pos & 63;

        ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64] =
            (unit & (((1ULL << writing_lower_bound) - 1) +
                     ((-1ULL) - ((-1ULL) >> (63 - writing_upper_bound))))) +
            ((all_encode &
              ((1ULL << (writing_upper_bound - writing_lower_bound + 1)) - 1))
             << writing_lower_bound);
    } else {
        uint64_t unit1 = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64];
        uint64_t unit2 = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64 + 1];
        int writing_lower_bound = start_bit_pos & 63;
        int writing_upper_bound = end_bit_pos & 63;
        uint64_t lower_part =
            all_encode & ((1LL << (64 - writing_lower_bound)) - 1);
        uint64_t higher_part = all_encode >> (64 - writing_lower_bound);
        ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64] =
            (unit1 & ((1LL << writing_lower_bound) - 1)) +
            (lower_part << writing_lower_bound);
        ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64 + 1] =
            ((unit2 >> (writing_upper_bound + 1)) << (writing_upper_bound + 1)) +
            higher_part;
    }

}


template <typename fp_t, int fp_len>
inline int SemiSortCuckooFilter<fp_t, fp_len>::high_bit(fp_t fp) {
    return (fp >> (fp_len - 4)) & ((1 << 4) - 1);
}

template <typename fp_t, int fp_len>
inline int SemiSortCuckooFilter<fp_t, fp_len>::low_bit(fp_t fp) {
    return fp & ((1 << (fp_len - 4)) - 1);
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::insert_to_bucket(fp_t* store, fp_t fp) {
    // if success return 0
    // if find collision : return 1 + position
    // if full : return 1 + 4

    if (store[this->m] == this->m)
        return 1 + 4;
    else {
        store[3] = fp; // sorted -- store[3] must be empty !
        return 0;
    }
}

template <typename fp_t, int fp_len>
bool SemiSortCuckooFilter<fp_t, fp_len>::insert(uint64_t ele) {
    return false;
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::lookup_in_bucket(int pos, fp_t fp) {
    // If lookup success return 1
    // If lookup fail and the bucket is full return 2
    // If lookup fail and the bucket is not full return 3

    fp_t store[8];
    get_bucket(pos, store);

    int isFull = 1;
    for (int i = 0; i < this->m; i++) {
        fp_t t = store[i];
        if (t == fp) return 1;
        isFull &= (t != 0);
    }
    return (isFull) ? 2 : 3;
}

template <typename fp_t, int fp_len>
bool SemiSortCuckooFilter<fp_t, fp_len>::lookup(uint64_t ele) {
    // If ele is positive, return true
    // negative -- return false

    fp_t fp = fingerprint(ele);
    int pos1 = this->position_hash(ele);

    int ok1 = lookup_in_bucket(pos1, fp);

    if (ok1 == 1) return true;
    // if (ok1 == 3) return false;

    int pos2 = alternate(pos1, fp);
    int ok2 = lookup_in_bucket(pos2, fp);

    return ok2 == 1;
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::del_in_bucket(int pos, fp_t fp) {
    fp_t store[8];
    get_bucket(pos, store);

    for (int i = 0; i < this->m; i++) {
        fp_t t = store[i];
        if (t == fp) {
            store[i] = 0;
            --this->filled_cell;
            set_bucket(pos, store);
            return 1;
        }
    }
    return 0;
}

template <typename fp_t, int fp_len>
bool SemiSortCuckooFilter<fp_t, fp_len>::del(uint64_t ele) {
    // If ele is positive, return true
    // negative -- return false

    fp_t fp = fingerprint(ele);
    int pos1 = this->position_hash(ele);

    int ok1 = del_in_bucket(pos1, fp);

    if (ok1 == 1) return true;
    // if (ok1 == 3) return false;

    int pos2 = alternate(pos1, fp);
    int ok2 = del_in_bucket(pos2, fp);

    return ok2 == 1;
}
template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_load_factor() {
    return filled_cell * 1.0 / this->n / this->m;
}


template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_bits_per_item() {
    return double(this->memory_consumption) * 8 / filled_cell;
}

template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_full_bucket_factor() {
    return full_bucket * 1.0 / this->n;
}
uint64_t _random_seed = 0x8091a2b3c4d5e6f7;
uint32_t hash_func3_32bit(uint32_t fp){
	uint64_t h = fp ^ _random_seed;
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}
template <typename fp_t, int fp_len>
int VacuumFilter<fp_t, fp_len>::alternate(int pos, fp_t fp)  // get alternate position
{
    uint32_t fp_hash = hash_func3_32bit(fp);
    int seg = this->len[fp & 3];
    return pos ^ (fp_hash & seg);
}

template <typename fp_t, int fp_len>
bool VacuumFilter<fp_t, fp_len>::insert(uint64_t ele) {
    // If insert success return true
    // If insert fail return false

    fp_t fp = this->fingerprint(ele);
    int cur1 = this->position_hash(ele);
    int cur2 = alternate(cur1, fp);

    fp_t store1[8];
    fp_t store2[8];

    this->get_bucket(cur1, store1);
    this->get_bucket(cur2, store2);

    if (store1[this->m] <= store2[this->m]) {
        if (this->insert_to_bucket(store1, fp) == 0) {
            this->filled_cell++;
            this->set_bucket(cur1, store1);
            return true;
        }
    } else {
        if (this->insert_to_bucket(store2, fp) == 0) {
            this->filled_cell++;
            this->set_bucket(cur2, store2);
            return true;
        }
    }

    // randomly choose one bucket's elements to kick
    int rk = rand() % this->m;

    // get those item
    int cur;
    fp_t* cur_store;

    if (rand() & 1)
        cur = cur1, cur_store = store1;
    else
        cur = cur2, cur_store = store2;

    fp_t tmp_fp = cur_store[rk];
    cur_store[rk] = fp;
    this->set_bucket(cur, cur_store);

    int alt = alternate(cur, tmp_fp);

    for (int i = 0; i < this->max_kick_steps; i++) {
        memset(store1, 0, sizeof(store1));
        this->get_bucket(alt, store1);
        if (store1[this->m] == this->m) {
            for (int j = 0; j < this->m; j++) {
                int nex = alternate(alt, store1[j]);
                this->get_bucket(nex, store2);
                if (store2[this->m] < this->m) {
                    store2[this->m - 1] = store1[j];
                    store1[j] = tmp_fp;
                    this->filled_cell++;
                    this->set_bucket(nex, store2);
                    this->set_bucket(alt, store1);
                    return true;
                }
            }

            rk = rand() % this->m;
            fp = store1[rk];
            store1[rk] = tmp_fp;
            this->set_bucket(alt, store1);

            tmp_fp = fp;
            alt = alternate(alt, tmp_fp);
        } else {
            store1[this->m - 1] = tmp_fp;
            this->filled_cell++;
            this->set_bucket(alt, store1);
            return true;
        }
    }
    return false;
}

template <typename fp_t, int fp_len>
VacuumFilter<fp_t, fp_len>::VacuumFilter(size_t nbits){
    this->init(nbits / (fp_len - 1), 4, 400);
    // mt19937 gen(1337);
    seed = rand(); //gen();
}
template <typename fp_t, int fp_len>
void VacuumFilter<fp_t, fp_len>::AddKeys(const vector<Bitwise> &keys){
    for(auto &key: keys)
        if(!this->insert(key.hash(seed)))
            throw std::runtime_error("Cuckoo Insert Failed");
}
template <typename fp_t, int fp_len>
void VacuumFilter<fp_t, fp_len>::AddKeys_len(const vector<Bitwise> &keys){
    for(auto &key: keys)
        if(!this->insert(key.hash_len(seed)))
            throw std::runtime_error("Cuckoo Insert Failed");
}
template <typename fp_t, int fp_len>
bool VacuumFilter<fp_t, fp_len>::Query(const Bitwise &key){
    return this->lookup(key.hash(seed));
}
template <typename fp_t, int fp_len>
bool VacuumFilter<fp_t, fp_len>::Query_len(const Bitwise &key){
    return this->lookup(key.hash_len(seed));
}
template <typename fp_t, int fp_len>
size_t VacuumFilter<fp_t, fp_len>::mem() const{
    return this->memory_consumption;
}

}


template class Rosetta<BloomFilter<false>, false>;
template class Rosetta<BloomFilter<true>, true>;
template class Rosetta<vacuum::VacuumFilter<uint16_t, 16>, false>;
template class Rosetta<vacuum::VacuumFilter<uint8_t, 8>, false>;
template class UnlayeredRosetta<vacuum::VacuumFilter<uint16_t, 16>, false>;
template class UnlayeredRosetta<vacuum::VacuumFilter<uint8_t, 8>, false>;

#ifdef USE_DTL
template class Rosetta<DtlBlockedBloomFilter>;
#endif
