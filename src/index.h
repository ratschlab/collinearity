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
#include "config.h"
#include "templated_tiered.h"

#ifdef NDEBUG
#define SANITY_CHECKS 0
#else
#define SANITY_CHECKS 1
#endif


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
    const u4 k, sigma, n_keys;
    const bool fwd_rev;
    const float presence_fraction;
    const int bandwidth;
    const u8 sort_blocksz;

    std::vector<std::string> headers;
    u4 max_occ = -1;
    parlay::sequence<u8> value_offsets;

    index_t();

public:
    index_t(config_t &config):
        k(config.k), sigma(config.sigma), fwd_rev(config.fwd_rev), sort_blocksz(config.sort_block_size),
        presence_fraction(config.presence_fraction), bandwidth(config.bandwidth), n_keys(1<<(config.k<<1)) {}
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
        if (seq.length() > 2 * k) {
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
        } else return {"*", true, 0, 0.0f};
    }
    /**
     * Build the index after adding all reference sequences
     */
    virtual void build() = 0;

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
    cqueue_t<u4> q_keys;
    cqueue_t<u4> q_values;
    heavyhitter_ht_t<u4> *hhs = nullptr;
    std::vector<u4> frag_offsets = {0};
    u4 frag_len, frag_ovlp_len;

    inline std::pair<cqueue_t<u4>::iterator_t, cqueue_t<u4>::iterator_t> get(u4 key) {
        return { q_values.it(value_offsets[key]), q_values.it(value_offsets[key+1]) };
    }

public:
    explicit j_index_t(config_t &config): index_t(config), frag_len(config.jc_frag_len), frag_ovlp_len(config.jc_frag_ovlp_len) {}
    void add(std::string &name, parlay::slice<char*, char*> seq) override;
    std::tuple<const char*, u4, float> search(parlay::slice<char*, char*> seq) override;
    void init_query_buffers() override;
    void build() override;
    void dump(const std::string& filename) override;
    void load(const std::string& filename) override;
};

class c_index_t : public index_t {
    cqueue_t<u4> q_keys;
    cqueue_t<u8> q_values;
    heavyhitter_ht_t<u8> *hhs = nullptr;

    inline std::pair<cqueue_t<u8>::iterator_t, cqueue_t<u8>::iterator_t> get(u4 key) {
        return { q_values.it(value_offsets[key]), q_values.it(value_offsets[key+1]) };
    }

public:
    explicit c_index_t(config_t &config) : index_t(config) {}
    void add(std::string &name, parlay::slice<char*, char*> seq) override;
    std::tuple<const char*, u4, float> search(parlay::slice<char*, char*> seq) override;
    void init_query_buffers() override;
    void build() override;
    void dump(const std::string& filename) override;
    void load(const std::string& filename) override;
};

#define N_SHARDS(n_shard_bits)                  (1<<n_shard_bits)
#define N_KEYS_PER_SHARD(n_keys, n_shard_bits)  ((n_keys) >> (n_shard_bits))
#define SHARD(x, n_shard_bits)                  ((x) & (N_SHARDS(n_shard_bits)-1))
#define SKEY(x, n_shard_bits)                   ((x) >> n_shard_bits)

using namespace Seq;
typedef Tiered<u8, LayerItr<LayerEnd,Layer<512,Layer<512,Layer<512>>>>> tvec_t;

class dindex_t {
    const int k, sigma, bandwidth, n_shard_bits, n_keys, n_shards, n_keys_per_shard;
    float presence_fraction;
    parlay::sequence<u4> keys;
    parlay::sequence<u8> values;
    heavyhitter_ht_t<u8> *hhs = nullptr;
    struct headers_t {
        std::vector<std::string> names;
        emhash8::HashMap<std::string, u8> name_map; /// header, id, length
        inline std::pair<u4, u4> get_id_offset(std::string &name) {
            if (name_map.contains(name)) {
                u8 val = name_map[name];
                u4 id = HIGH32(val), offset = LOW32(val);
                return {id, offset};
            } else {
                u4 id = names.size();
                name_map[name] = MAKE64(id, 0);
                names.push_back(name);
                return {id, 0};
            }
        }
        inline std::string& get_name(u4 id) { return names[id]; }
    };
    headers_t headers;
    struct shard_t {
        tvec_t values;
        std::vector<u8> offsets;
        explicit shard_t(size_t n_keys) {
            offsets.resize(n_keys);
            p_fill(offsets.data(), offsets.size(), 0ul);
        }
    };
    class iterator_t {
    private:
        tvec_t &v;
        size_t i;
    public:
        using value_type = u8;
        using difference_type = size_t;
        using pointer = u8*;
        using reference = u8&;
        using iterator_category = std::forward_iterator_tag;

        // Constructor
        explicit iterator_t(tvec_t &v, size_t i=0) : v(v), i(i) {}

        // Dereference operator
        u8 operator*() { return v[i]; }

        // Arrow operator
        const u8* operator->() { return &(v[i]); }

        // Prefix increment
        iterator_t& operator++() {
            ++i;
            return *this;
        }

        // Postfix increment
        iterator_t operator++(int) {
            iterator_t temp = *this;
            ++(*this);
            return temp;
        }

        // Comparison operators
        bool operator==(const iterator_t& other) const { return i == other.i; }
        bool operator!=(const iterator_t& other) const { return i != other.i; }
        bool operator<(const iterator_t& other) const { return i < other.i; }
        bool operator<=(const iterator_t& other) const { return i <= other.i; }
        bool operator>(const iterator_t& other) const { return i > other.i; }
        bool operator>=(const iterator_t& other) const { return i >= other.i; }
    };
    std::vector<shard_t> shards;
    void put_in_shard(shard_t &shard, parlay::sequence<u4> &keys, parlay::sequence<u8> &values, parlay::sequence<u8> &indices);
    std::pair<iterator_t,iterator_t> get(u4 key);

public:

    explicit dindex_t(config_t &config): k(config.k), sigma(config.sigma),
                                         presence_fraction(config.presence_fraction), bandwidth(config.bandwidth),
                                         n_keys(1<<(config.k<<1)), n_shard_bits(config.n_shard_bits),
                                         n_shards(N_SHARDS(n_shard_bits)), n_keys_per_shard(N_KEYS_PER_SHARD(n_keys, n_shard_bits)) {
        for (int i = 0; i < n_shards; ++i)
            shards.emplace_back(n_keys_per_shard+1);
        hhs = new heavyhitter_ht_t<u8>[parlay::num_workers()];
    }

    void add(std::string &name, parlay::slice<char*, char*> seq);
    std::tuple<const char*, u4, float> search(parlay::slice<char*, char*> seq);
    void merge();
    std::tuple<const char*, bool, u4, float> search(std::string &seq);
};

#endif //COLLINEARITY_INDEX_H
