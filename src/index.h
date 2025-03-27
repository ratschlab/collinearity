//
// Created by Sayan Goswami on 09.12.2024.
//

#ifndef COLLINEARITY_INDEX_H
#define COLLINEARITY_INDEX_H

#include "prelude.h"
#include "parlay_utils.h"
#include "cqueue.h"
#include "cqutils.h"
#include "hash_table8.hpp"
#include "utils.h"

#ifdef NDEBUG
#define SANITY_CHECKS 0
#else
#define SANITY_CHECKS 1
#endif

#define N_SHARD_BITS 7
#define N_SHARDS (1<<N_SHARD_BITS)
#define SHARD(x) ((x) & (N_SHARDS-1))
#define SKEY(x) ((x) >> N_SHARD_BITS)

// the minimum fraction of kmers required for majority vote
const float presence_fraction = .1f;

// this is the bandwidth, I will later set this from the config
const int bandwidth = 15;

/**
 * A shard which contains keys and a sequence of values for each key
 * @tparam K key type
 * @tparam V value type
 */
template <typename K, typename V>
struct shard_t {
    emhash8::HashMap<K,u8> tuples;
    parlay::sequence<V> values;
};

/**
 * A simple frequency counter
 * @tparam T key type
 */
template <typename T>
struct heavyhitter_ht_t {
    emhash8::HashMap<T,u4> counts;
    T top_key = -1;
    u4 top_count = 0;
    void insert(const T key) {
        u4 count = counts[key]++;
        if (count > top_count)
            top_count = count, top_key = key;
    }
    void reset() { counts.clear(), top_key=-1, top_count = 0; }
};

/**
 * An interface for an index
 */
class index_t {
protected:
    const u4 k, sigma;
    const bool fwd_rev;

    std::vector<std::string> headers;
    u4 max_occ = -1;

    index_t();

public:
    index_t(u4 k, u4 sigma, bool fwd_rev): k(k), sigma(sigma), fwd_rev(fwd_rev) {}
    virtual ~index_t() {}

    /**
     * Initialize buffers for query client
     */
    virtual void init_query_buffers() = 0;

    /**
     * Add a sequence to the index
     * @param name reference header
     * @param seq a parlay slice view of a reference sequence
     */
    virtual void add(std::string &name, parlay::slice<char*, char*> seq) = 0;

    /**
     * Search for a sequence in the index
     * @param seq a parlay slice view of a query sequence
     * @return a tuple consisting of [
     * (1) the reference header or * if no match was found,
     * (2) position of query in reference or 2^32 if no match was found, and
     * (3) fraction of k-mers supporting the match
     * ]
     */
    virtual std::tuple<const char*, u4, float> search(parlay::slice<char*, char*> seq) = 0;

    /**
     * Add a sequence to the index
     * @param name reference header
     * @param seq reference sequence
     */
    void add(std::string &name, std::string &seq) {
        if (fwd_rev) {
            auto s_name = name + "+";
            add(s_name, parlay::make_slice(seq.data(), seq.data() + seq.size()));
            s_name = name + "-";
            auto rev = revcmp(seq);
            add(s_name, parlay::make_slice(rev.begin(), rev.end()));
        } else {
            add(name, parlay::make_slice(seq.data(), seq.data() + seq.size()));
        }
    }

    /**
     * Search for a sequence in the index
     * @param seq query sequence
     * @return a tuple consisting of [
     * (1) the reference header or * if no match was found,
     * (2) true if the query aligned with forward strand and false otherwise,
     * (3) position of query in reference or 2^32 if no match was found, and
     * (4) fraction of k-mers supporting the match
     * ]
     */
    std::tuple<const char*, bool, u4, float> search(std::string &seq) {
        const auto [header1, pos1, support1] =
                search(parlay::make_slice(seq.data(), seq.data() + seq.size()));
        if (fwd_rev) return std::make_tuple(header1, true, pos1, support1);
        else {
            auto slice = parlay::make_slice(seq.rbegin(), seq.rend());
            auto rc = parlay::map(slice, [](char c) {
                return "TGAC"[(c >> 1) & 3];
            });
            const auto [header2, pos2, support2] =
                    search(parlay::make_slice(rc.begin(), rc.end()));
            if (support1 >= support2)
                return std::make_tuple(header1, true, pos1, support1);
            else
                return std::make_tuple(header2, false, pos2, support2);
        }
    }
    /**
     * Build the index after adding all reference sequences
     */
    virtual void build() = 0;
//    std::tuple<const char*, bool, u4, float> get(std::string &seq) = 0;
//    {
//        bool qry_fwd = true;
//
//        const auto i = parlay::worker_id();
//        auto &hh = hhs[i];
//
//        parlay::sequence<u4> keys = create_kmers_1t(seq, k, sigma, encode_dna);
//        hh.reset();
//
//        for (u4 j = 0; j < keys.size(); ++j) {
//            const auto &[vbegin, vend] = get(keys[j]);
//            for (auto v = vbegin; v != vend; ++v) update_count(hh, j, *v);
//        }
//
//        auto top_count = hh.top_count;
//        auto top_key = hh.top_key;
//
//        if (!fwd_rev) { // only one strand is indexed. must query both strands
//            auto slice = parlay::make_slice(seq.rbegin(), seq.rend());
//            auto rc = parlay::delayed_map(slice, [](char c) {
//                return "TGAC"[(c >> 1) & 3];
//            });
//            keys = create_kmers_1t(rc, k, sigma, encode_dna);
//            hh.reset();
//
//
//            for (u4 j = 0; j < keys.size(); ++j) {
//                const auto &[vbegin, vend] = get(keys[j]);
//                for (auto v = vbegin; v != vend; ++v) update_count(hh, j, *v);
//            }
//
//            if (hh.top_count > top_count) {
//                top_count = hh.top_count;
//                top_key = hh.top_key;
//                qry_fwd = false;
//            }
//        }
//
//        // presence_fraction is the frac. of kmers that support a match (or intercept)
//        float presence = (top_count * 1.0) / keys.size();
//        if (presence >= presence_fraction) return get_alignment(top_key, presence, qry_fwd);
//        else return std::make_tuple("*", true, -1, 0.0f);
//    }

    /**
     * Dump index into file
     * @param filename name (path) of dump file
     */
    virtual void dump(const std::string& filename) = 0;

    /**
     * Load index from file
     * @param filename name (path) of dump file
     */
    virtual void load(const std::string& filename) = 0;

};

/**
 * index to detect jaccard similarity
 */
class j_index_t : public index_t {
    shard_t<u4,u4> shards[N_SHARDS];
    cqueue_t<u4> q_keys;
    cqueue_t<u4> q_values;
    heavyhitter_ht_t<u4> *hhs = nullptr;
    std::vector<u4> frag_offsets = {0};
    u4 frag_len, frag_ovlp_len;

public:
    j_index_t(u4 k, u4 sigma, bool fwd_rev, u4 frag_len=180, u4 frag_ovlp_len=120):
        index_t(k, sigma, fwd_rev), frag_len(frag_len), frag_ovlp_len(frag_ovlp_len) {}
    void add(std::string &name, parlay::slice<char*, char*> seq) override;
    std::tuple<const char*, u4, float> search(parlay::slice<char*, char*> seq) override;
    void init_query_buffers() override;
    void build() override;
    void dump(const std::string& filename) override;
    void load(const std::string& filename) override;
};

class c_index_t : public index_t {
    shard_t<u4,u8> shards[N_SHARDS];
    cqueue_t<u4> q_keys;
    cqueue_t<u8> q_values;
    heavyhitter_ht_t<u8> *hhs = nullptr;
public:
    c_index_t(u4 k, u4 sigma, bool fwd_rev) : index_t(k, sigma, fwd_rev) {}
    void add(std::string &name, parlay::slice<char*, char*> seq) override;
    std::tuple<const char*, u4, float> search(parlay::slice<char*, char*> seq) override;
    void init_query_buffers() override;
    void build() override;
    void dump(const std::string& filename) override;
    void load(const std::string& filename) override;
};

template <typename K, typename V>
static u4 consolidate(cqueue_t<K> &q_keys, cqueue_t<V> &q_values, shard_t<K, V> *shards);

//struct sindex_t {
//    shard_t<u4,u8> shards[N_SHARDS];
//    cqueue_t<u4> q_keys;
//    cqueue_t<u8> q_values;
//    std::vector<std::string> headers;
//    u4 max_allowed_occ;
//    heavyhitter_ht_t<u8> *hhs = nullptr;
//    const bool minimize = true;
//
//    void add(std::string &name, parlay::sequence<u4> &keys, parlay::sequence<u8> &addresses) {
//        headers.push_back(name);
//        q_keys.push_back(keys.data(), keys.size());
//        q_values.push_back(addresses.data(), addresses.size());
//    }
//
//    void build() {
//        max_allowed_occ = consolidate(q_keys, q_values, shards);
//    }
//
//    void init_query_buffers() {hhs = new heavyhitter_ht_t<u8>[parlay::num_workers()];}
//
//    std::tuple<u4, u4, float> search(const u4 *keys, const size_t n_keys) {
//        const auto i = parlay::worker_id();
//        auto &hh = hhs[i];
//        hh.reset();
//
//        for (u4 j = 0; j < n_keys; ++j) {
//            const auto &[vbegin, vend] = get(keys[j]);
//            for (auto v = vbegin; v != vend; ++v) {
//                u8 ref_id = get_id_from(*v);
//                u8 ref_pos = get_pos_from(*v);
//                u8 intercept = (ref_pos > j) ? (ref_pos - j) : 0;
//                intercept /= bandwidth;
//                u8 key = make_key_from(ref_id, intercept);
//                hh.insert(key);
//                if (intercept >= bandwidth) {
//                    intercept -= bandwidth;
//                    key = make_key_from(ref_id, intercept);
//                    hh.insert(key);
//                }
//            }
//        }
//
//        // in this case, presence_fraction is the frac. of kmers that support an intercept
//        float presence = (hh.top_count * 1.0) / n_keys;
//        if (presence >= presence_fraction)
//            return std::make_tuple(get_id_from(hh.top_key), get_pos_from(hh.top_key) * bandwidth, presence);
//        else return std::make_tuple(-1, -1, 0.0f);
//    }
//
//    std::pair<u8*, u8*> get(u4 key) {
//        u4 shard_id = SHARD(key);
//        auto &shard = shards[shard_id];
//        if (shard.tuples.contains(SKEY(key))) {
//            auto val = shard.tuples[SKEY(key)];
//            auto offset = val >> 32, count = val & 0xffffffff;
////            if (count >= max_allowed_occ) return {nullptr, nullptr};
//            return {shard.values.data() + offset, shard.values.data() + offset + count};
//        } else return {nullptr, nullptr};
//    }
//
//    void dump(const std::string& basename) {
//        log_info("Dumping index.");
//        std::string filename = basename + ".cidx";
//        FILE *fp = fopen(filename.c_str(), "wb");
//        if (!fp) log_error("Could not open %s because %s", filename.c_str(), strerror(errno));
//
//        size_t n = headers.size();
//        log_info("Dumping %zd refnames..", n);
//        fwrite(&n, sizeof(n), 1, fp);
//        std::vector<u2> refnamelens(n);
//        for (int i = 0; i < n; ++i) refnamelens[i] = headers[i].size();
//        fwrite(refnamelens.data(), sizeof(u2), n, fp);
//        for (const auto& refname : headers) {
//            n = refname.size();
//            fwrite(refname.c_str(), n, 1, fp);
//        }
//
//        log_info("Collecting key-value pairs from shards..");
//        parlay::sequence<u4> n_unique_keys_per_shard(N_SHARDS);
//        parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
//            n_unique_keys_per_shard[i] = shards[i].tuples.size();
//        });
//        u4 total_uniq_keys = parlay::reduce(n_unique_keys_per_shard);
//        parlay::sequence<u4> keys(total_uniq_keys);
//        parlay::sequence<u8> values(total_uniq_keys);
//        auto [offset, total] = parlay::scan(n_unique_keys_per_shard);
//        expect(total == total_uniq_keys);
//        offset.push_back(total_uniq_keys);
//
//        parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
//            auto my_off = offset[i];
//            int j = 0;
//            for (auto &[key, value]: shards[i].tuples) {
//                keys[my_off + j] = key;
//                values[my_off + j] = value;
//                j++;
//            }
//            expect(j == offset[i+1] - offset[i]);
//        });
//
//        log_info("Dumping keys..");
//        expect(offset.size() == N_SHARDS + 1);
//        fwrite(&total_uniq_keys, sizeof(total_uniq_keys), 1, fp);
//        fwrite(offset.data(), sizeof(offset[0]), offset.size(), fp);
//        fwrite(keys.data(), sizeof(keys[0]), keys.size(), fp);
//        fwrite(values.data(), sizeof(values[0]), values.size(), fp);
//
//        log_info("Dumping values..");
//        for (auto & shard : shards) {
//            u4 n_values = shard.values.size();
//            fwrite(&n_values, sizeof(n_values), 1, fp);
//            fwrite(shard.values.data(), sizeof(shard.values[0]), shard.values.size(), fp);
//        }
//
//        fwrite(&max_allowed_occ, sizeof(max_allowed_occ), 1, fp);
//
//        fclose(fp);
//        log_info("Done.");
//    }
//
//    void load(const std::string& basename) {
//        log_info("Loading index.");
//        std::string filename = basename + ".cidx";
//        FILE *fp = fopen(filename.c_str(), "rb");
//        if (!fp) log_error("Could not open %s because %s", filename.c_str(), strerror(errno));
//
//        size_t n;
//        expect(fread(&n, sizeof(n), 1, fp) == 1);
//        log_info("Loading %zd refnames..", n);
//        std::vector<u2> refnamelens(n);
//        expect(fread(refnamelens.data(), sizeof(u2), n, fp) == n);
//
//        char buffer[4096];
//        for (int i = 0; i < n; ++i) {
//            size_t l = refnamelens[i];
//            expect(fread(buffer, l, 1, fp) == 1);
//            headers.emplace_back(buffer, l);
//        }
//
//        log_info("loading keys..");
//        u4 total_uniq_keys;
//        parlay::sequence<u4> offset(N_SHARDS+1);
//        expect(fread(&total_uniq_keys, sizeof(total_uniq_keys), 1, fp) == 1);
//        expect(fread(offset.data(), sizeof(offset[0]), offset.size(), fp) == offset.size());
//
//        parlay::sequence<u4> keys(total_uniq_keys);
//        parlay::sequence<u8> values(total_uniq_keys);
//        expect(fread(keys.data(), sizeof(keys[0]), keys.size(), fp) == keys.size());
//        expect(fread(values.data(), sizeof(values[0]), values.size(), fp) == values.size());
//
//        log_info("adding keys into hash table");
//        parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
//            auto my_off = offset[i], count = offset[i+1] - offset[i];
//            for (int j = 0; j < count; ++j)
//                shards[i].tuples[keys[my_off + j]] = values[my_off + j];
//        });
//
//        log_info("loading values..");
//        for (auto & shard : shards) {
//            u4 n_values = shard.values.size();
//            expect(fread(&n_values, sizeof(n_values), 1, fp) == 1);
//            shard.values.resize(n_values);
//            expect(fread(shard.values.data(), sizeof(shard.values[0]), shard.values.size(), fp) == shard.values.size());
//        }
//
//        expect(fread(&max_allowed_occ, sizeof(max_allowed_occ), 1, fp) == 1);
//
//        fclose(fp);
//        log_info("Done.");
//    }
//};

#endif //COLLINEARITY_INDEX_H
