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
    std::vector<compressed_array_t<u8>> values;

    index_t(CQueue<u4> &qkeys, CQueue<u4> &qcounts, CQueue<u8> &qvalues) {
        values.resize(NKEYS);
        size_t nk = qkeys.size();
        expect(nk == qcounts.size());
        std::vector<u4> keys(nk), counts(nk);
        qkeys.pop_front(keys.data(), nk);
        qcounts.pop_front(counts.data(), nk);
        verify(parlay::reduce(counts) == qvalues.size());

        info("Consolidating values for each unique key.");
        size_t max_key_count = *parlay::max_element(counts);
        std::vector<u8> val_buf(max_key_count);
        for (int i = 0; i < nk; ++i) {
            auto key = keys[i];
            auto count = counts[i];
            size_t nv = qvalues.pop_front(val_buf.data(), count);
            verify(nv == count);
            values[key].compress(val_buf.data(), nv, true);
        }
        info("Done");
    }

    const std::vector<u8> &get(u4 key) {
        return values[key].v;
    }
};

#endif //COLLINEARITY_INDEX_H
