//
// Created by Sayan Goswami on 29.11.2024.
//

#include "tests.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

int main(int argc, char *argv[]) {
//    test_cqueue();
//    test_cqutils();
    log_info("Number of thread ids assigned - %d", parlay::num_thread_ids());
    log_info("Number of workers - %zd", parlay::num_workers());
    parlay::for_each(parlay::iota(32), [](size_t i){
        printf("%zd - %zd\n", i, parlay::worker_id());
    });
    log_info("Number of thread ids assigned - %d", parlay::num_thread_ids());
//    test_compressed_array();
//    test_cqueue();
//    test_cqutils();

    // Input: a sorted sequence
//    parlay::sequence<int> sorted_seq = {0,1,1,2,2,2,2,3,3,3,4,4,4,5,6,6,7,8,8,8,9};
//
//    auto kv = parlay::map(sorted_seq, [&](size_t i){
//        return std::make_pair(sorted_seq[i], 1);
//    });
//    auto rle = parlay::group_by_key_ordered(kv);
//
//    auto histogram = parlay::histogram_by_key(sorted_seq);
//    parlay::sort_inplace(histogram);
//    for (auto &[c, n] : histogram)
//        printf("%d - %zd\n", c, n);
}