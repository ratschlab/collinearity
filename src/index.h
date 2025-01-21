//
// Created by Sayan Goswami on 09.12.2024.
//

#ifndef COLLINEARITY_INDEX_H
#define COLLINEARITY_INDEX_H

#include "prelude.h"
#include "compressed_array.h"
#include "parlay_utils.h"
//#include "parlay_hash/unordered_map.h"
#include "cqueue.h"
#include "tsl/hopscotch_map.h"
#include "hash_table8.hpp"

#define SANITY_CHECKS 1

#define NKEYS (1<<(KMER_LENGTH<<1))

struct index_t {
    std::vector<u4> counts;
    std::vector<u8> offsets, values;
    index_t() = default;
    index_t(cqueue_t<u4> &qkeys, cqueue_t<u4> &qcounts, cqueue_t<u8> &qvalues) {
        values.resize(NKEYS);
        size_t nk = qkeys.size();
        size_t nv = qvalues.size();
        expect(nk == qcounts.size());
        std::vector<u4> keys(nk), n_values(nk);
        qkeys.pop_front(keys.data(), nk);
        qcounts.pop_front(n_values.data(), nk);
        if (SANITY_CHECKS) {
            size_t total = 0;
            for (auto c: n_values) total += c;
            expect(total == nv);
        }

        info("Allocating memory for coordinates..");
        values.resize(nv);

        info("Consolidating coordinates..");
        qvalues.pop_front(values.data(), values.size());
        counts.resize(NKEYS);
        offsets.resize(NKEYS);
        std::fill(counts.begin(), counts.end(), 0);
        std::fill(offsets.begin(), offsets.end(), 0);
        parlay::for_each(parlay::iota(nk), [&](size_t i){
            auto key = keys[i];
            counts[key] = n_values[i];
            offsets[key] = n_values[i];
        });
        size_t c = parlay::scan_inplace(offsets);
        verify(c == values.size());
        info("Done");
    }

    std::pair<u8*, u8*> get(u4 key) {
        return std::make_pair(values.data()+offsets[key], values.data()+offsets[key]+counts[key]);
    }
};

struct [[maybe_unused]] heavyhitter_ht_t {
    emhash8::HashMap<u8,u4> counts;
    u8 top_key = -1;
    u4 top_count = 0;
    void insert(const u8 key) {
        u4 count = counts[key]++;
        if (count > top_count)
            top_count = count, top_key = key;
    }
    void reset() { counts.clear(), top_key=-1, top_count = 0; }
};

/**
 * dump index into file
 * @param idx index
 * @param refnames reference headers
 * @param basename the filepath where the index will be dumped will be "{basename}.cidx"
 */
static void dump_index(index_t *idx, std::vector<std::string> &refnames, const char *basename) {
    info("Dumping index.");
    std::string filename = std::string(basename) + std::string(".cidx");
    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp) error("Could not open %s because %s", filename.c_str(), strerror(errno));

    size_t n = refnames.size();
    info("Dumping %zd refnames..", n);
    fwrite(&n, sizeof(n), 1, fp);
    std::vector<u2> refnamelens(n);
    for (int i = 0; i < n; ++i) refnamelens[i] = refnames[i].size();
    fwrite(refnamelens.data(), sizeof(u2), n, fp);
    for (const auto& refname : refnames) {
        n = refname.size();
        fwrite(refname.c_str(), n, 1, fp);
    }

    n = idx->counts.size();
    info("Dumping %zd counts and offsets..", n);
    fwrite(&n, sizeof(n), 1, fp);
    fwrite(idx->counts.data(), sizeof(u4), n, fp);
    fwrite(idx->offsets.data(), sizeof(u8), n, fp);

    n = idx->values.size();
    info("Dumping %zd values..", n);
    fwrite(&n, sizeof(n), 1, fp);
    fwrite(idx->values.data(), sizeof(u8), n, fp);
    fclose(fp);
    info("Done.");
}

/**
 * load index from file
 * @param filename path to the .cidx file where the index was dumped
 * @return a pair of index and headers
 */
static std::pair<index_t*, std::vector<std::string>> load_index(char *filename) {
    info("Loading index.");
    FILE *fp = fopen(filename, "rb");
    if (!fp) error("Could not open %s because %s", filename, strerror(errno));

    size_t n;
    fread(&n, sizeof(n), 1, fp);
    info("Loading %zd refnames..", n);
    std::vector<u2> refnamelens(n);
    fread(refnamelens.data(), sizeof(u2), n, fp);

    std::vector<std::string> refnames;
    char buffer[4096];
    for (int i = 0; i < n; ++i) {
        size_t l = refnamelens[i];
        fread(buffer, l, 1, fp);
        refnames.emplace_back(buffer, l);
    }

    auto *idx = new index_t;
    fread(&n, sizeof(n), 1, fp);
    info("Loading %zd counts and offsets..", n);
    idx->counts.resize(n);
    idx->offsets.resize(n);
    fread(idx->counts.data(), sizeof(u4), n, fp);
    fread(idx->offsets.data(), sizeof(u8), n, fp);

    fread(&n, sizeof(n), 1, fp);
    info("Loading %zd values..", n);
    idx->values.resize(n);
    fread(idx->values.data(), sizeof(u8), n, fp);

    fclose(fp);
    return std::make_pair(idx, std::move(refnames));
}

#endif //COLLINEARITY_INDEX_H
