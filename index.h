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
// K is hardcoded here. change later
#define K 10
#define NKEYS (1<<(K<<1))

struct index_t {
    std::vector<u4> counts, offsets;
    std::vector<u8> values;

    index_t(CQueue<u4> &qkeys, CQueue<u4> &qcounts, CQueue<u8> &qvalues) {
        values.resize(NKEYS);
        size_t nk = qkeys.size();
        expect(nk == qcounts.size());
        std::vector<u4> keys(nk), n_values(nk);
        qkeys.pop_front(keys.data(), nk);
        qcounts.pop_front(n_values.data(), nk);
        verify(parlay::reduce(n_values) == qvalues.size());

        info("Allocating memory for coordinates..");
        values.resize(qvalues.size());

        info("Consolidating coordinates..");
        qvalues.pop_front(values.data(), values.size());
        counts.resize(NKEYS);
        offsets.resize(NKEYS);
        std::fill(counts.begin(), counts.end(), 0);
        std::fill(offsets.begin(), offsets.end(), 0);
        parlay::for_each(parlay::iota(nk), [&](size_t i){
            auto key = keys[i];
            counts[key] = n_values[key];
            offsets[key] = n_values[key];
        });
        parlay::scan_inplace(offsets);
        info("Done");
    }

    std::pair<u8*, u8*> get(u4 key) {
        return std::make_pair(values.data()+offsets[key], values.data()+offsets[key]+counts[key]);
    }
};

#endif //COLLINEARITY_INDEX_H
