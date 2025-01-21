//
// Created by Sayan Goswami on 29.11.2024.
//

#include "tests.h"
#include "../src/cqutils.h"

void test_cqutils() {
    const size_t N = (4 MiB) + 123456;
    auto keys = generate_random<u4>(N<<4, N);
    auto values = generate_random<u8>(N<<4, N);

    cqueue_t<u4> cqkeys(1 MiB);
    cqkeys.push_back((const u4*)keys.data(), N);

    cqueue_t<u8> cqvalues(1 MiB);
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
    check(parlay::is_sorted(res_keys), "Result keys are sorted");

    cqkeys.push_back(res_keys.data(), N);
    cqueue_t<u4> unique_keys(1 MiB);
    cqueue_t<u4> counts(1 MiB);
    cq_count_unique(cqkeys, 1 MiB, buf, unique_keys, counts);
    size_t Nu = unique_keys.size();
    info("Found %zd unique keys", Nu);
    check(counts.size() == Nu, "Each key has a count");
    parlay::sequence<u4> unique_keys_sorted(Nu), unique_key_counts(Nu);
    unique_keys.pop_front(unique_keys_sorted.data(), Nu);
    counts.pop_front(unique_key_counts.data(), Nu);
    check(parlay::is_sorted(unique_keys_sorted), "Unique keys are sorted");

    auto histogram = parlay::histogram_by_key(keys);
    parlay::sort_inplace(histogram);
    check(histogram.size() == Nu, "Correct count of keys");

    auto same = parlay::tabulate(Nu, [&](size_t i){
        return histogram[i].first == unique_keys_sorted[i] && histogram[i].second == unique_key_counts[i];
    });
    check(parlay::all_of((same), [&](auto &v){return v;}), "Correct keys and counts");

    unique_keys.clear();
    counts.clear();

}