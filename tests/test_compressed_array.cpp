//
// Created by Sayan Goswami on 29.11.2024.
//

#include "tests.h"
#include "../compressed_array.h"
#include <parlay/primitives.h>
#include <parlay/parallel.h>

void test_compressed_array() {
    const u4 N = 64<<20;
    const u4 I = 1<<10;

    // generate random numbers
    parlay::random_generator gen;
    std::uniform_int_distribution<u4> dis(0, N>>5), idis(0, N);
    auto input = parlay::tabulate(N, [&](size_t i) {
        auto r = gen[i];
        return dis(r);
    });
    auto rand_idx = parlay::tabulate(I, [&](size_t i) {
        auto r = gen[i];
        return idis(r);
    });

    // create compressed array
    compressed_array_t<u4, false> carr(input.data(), input.size());
    info("Compression ratio = %.2f%%", carr.nc * 100.0 / (carr.n * sizeof(u4)));

    // decompress and check for equality
    parlay::sequence<u4> buffer(N);
    carr.decompress(buffer.data());
    test(parlay::equal(input, buffer), "Compress/Decompress unsorted");

    // test selection unsorted for a random sample
    parlay::sequence<u4> diff(I);
    parlay::for_each(diff, [&](size_t i){
        diff[i] = input[rand_idx[i]] - carr[rand_idx[i]];
    });
    test(parlay::reduce(diff)==0, "Select unsorted");

    // sort input
    parlay::integer_sort_inplace(input);
    expect(parlay::is_sorted(input));

    // create compressed array
    compressed_array_t<u4, true> carr1(input.data(), input.size());
    info("Compression ratio = %.2f%%", carr1.nc * 100.0 / (carr1.n * sizeof(u4)));

    // decompress and check for equality
    carr1.decompress(buffer.data());
    test(parlay::equal(input, buffer), "Compress/Decompress sorted");

    // test selection sorted for a random sample
    parlay::for_each(diff, [&](size_t i){
        diff[i] = input[rand_idx[i]] - carr1[rand_idx[i]];
    });
    test(parlay::reduce(diff)==0, "Select sorted");
}
