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


#endif //COLLINEARITY_INDEX_H
