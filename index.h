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

#define SANITY_CHECKS 1

#define NKEYS (1<<(KMER_LENGTH<<1))

struct index_t {
    std::vector<u4> counts;
    std::vector<u8> offsets, values;
    index_t() = default;
    index_t(CQueue<u4> &qkeys, CQueue<u4> &qcounts, CQueue<u8> &qvalues) {
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

#endif //COLLINEARITY_INDEX_H
