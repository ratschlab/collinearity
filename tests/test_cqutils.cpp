//
// Created by Sayan Goswami on 29.11.2024.
//

#include "tests.h"
#include "../cqutils.h"

void test_cqutils() {
    const size_t N = (4 MiB) + 123456;
    auto keys = generate_random<u4>(N<<5, N);
    auto values = generate_random<u8>(N<<5, N);

    CQueue<u4> cqkeys(1 MiB);
    cqkeys.push_back((const u4*)keys.data(), N);

    CQueue<u8> cqvalues(1 MiB);
    cqvalues.push_back((const u8*)values.data(), N);

    auto buf = malloc(12 MiB);
    cq_sort_by_key(cqkeys, cqvalues, 1 MiB, buf);
    info("Sorted queue.");

    sort_by_key(keys.data(), values.data(), N);
    info("Sorted parlay");

    parlay::sequence<u4> res_keys(N);
    parlay::sequence<u8> res_values(N);
    cqkeys.pop_front(res_keys.data(), N);
    cqvalues.pop_front(res_values.data(), N);

    check(parlay::equal(keys, res_keys), "Keys match");
    check(parlay::equal(values, res_values), "Values match");


}