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
//            parlay::sort_inplace(val_buf);
            values[key].compress(val_buf.data(), nv, true);
        }
        info("Done");
    }

    const std::vector<u8> &get(u4 key) {
        return values[key].v;
    }

//    void batch_get(parlay::sequence<u4> &qry_keys, idx_slice_t &slice) {
//        slice.reset();
//        auto histogram = parlay::histogram_by_key(qry_keys);
//        parlay::sort_inplace(histogram);
////        size_t num_unique_qry_keys = histogram.size();
////        slice.val_counts.resize(num_unique_qry_keys);
////        slice.val_offsets.resize(num_unique_qry_keys);
//
////        parlay::for_each(parlay::iota(histogram.size()), [&](size_t i){
//        for (u8 i = 0; i < histogram.size(); ++i) {
//            auto key = histogram[i].first;
//            slice.val_counts[key] = values[key].n;
//            slice.val_offsets[key] = slice.val_counts[key];
//
////            auto val_p = kv.find(key);
////            if (val_p != kv.end()) {
////                auto idx = val_p.value();
////                slice.val_counts[i] = values[idx].n;
////                slice.val_offsets[i] = slice.val_counts[i];
////            } else {
////                slice.val_counts[i] = 0;
////                slice.val_offsets[i] = 0;
////            }
//////            slice.kv.Insert(key, i);
////            slice.kv[key] = i;
//        }
////        );
//        slice.values.resize(parlay::scan_inplace(slice.val_offsets));
//        expect(slice.values.size() == parlay::reduce(slice.val_counts));
//
////        parlay::for_each(parlay::iota(histogram.size()), [&](size_t i){
//        const auto valptr = slice.values.data();
//        for (u8 i = 0; i < histogram.size(); ++i) {
//            auto key = histogram[i].first;
//            if (slice.val_counts[key]) {
//
////                auto idx = kv[key];
//
//                size_t n = values[key].decompress(valptr + slice.val_offsets[key]);
//                verify(n == slice.val_counts[key]);
//            }
//        }
//        );
//    }
};

#endif //COLLINEARITY_INDEX_H
