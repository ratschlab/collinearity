//
// Created by Sayan Goswami on 29.11.2024.
//

#include "tests.h"
#include "../compressed_array.h"
#include <parlay/primitives.h>
#include <parlay/parallel.h>

void test_compressed_array() {
    const u4 N = 1<<20;
    const u4 I = 1<<6;

    // generate random numbers
    auto input = generate_random<u4>(N>>5, N);
    auto rand_idx = generate_random<u4>(N, I);

    // create compressed array
    compressed_array_t<u4> carr(input.data(), input.size());
//    info("Compression ratio = %.2f%%", carr.nc * 100.0 / (carr.n * sizeof(u4)));

    // decompress and check for equality
    parlay::sequence<u4> buffer(N);
    carr.decompress(buffer.data());
    check(parlay::equal(input, buffer), "Compress/Decompress unsorted");

    // test selection unsorted for a random sample
    parlay::sequence<u4> diff(I);
    parlay::for_each(diff, [&](size_t i){
        diff[i] = input[rand_idx[i]] - carr[rand_idx[i]];
    });
    check(parlay::reduce(diff)==0, "Select unsorted");

    // sort input
    parlay::integer_sort_inplace(input);
    expect(parlay::is_sorted(input));

    // create compressed array
    compressed_array_t<u4> carr1(input.data(), input.size(), true);
//    info("Compression ratio = %.2f%%", carr1.nc * 100.0 / (carr1.n * sizeof(u4)));

    // decompress and check for equality
    carr1.decompress(buffer.data());
    check(parlay::equal(input, buffer), "Compress/Decompress sorted");

    // test selection sorted for a random sample
    parlay::for_each(diff, [&](size_t i){
        diff[i] = input[rand_idx[i]] - carr1[rand_idx[i]];
    });
    check(parlay::reduce(diff)==0, "Select sorted");
}
