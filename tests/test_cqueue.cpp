//
// Created by Sayan Goswami on 29.11.2024.
//

#include "tests.h"
#include "../cqueue.h"

void test_cqueue() {
    const size_t N = (16 MiB) + 123456;
    auto input = generate_random<u4>(N<<5, N);
    parlay::sequence<u4> buffer(N);

    CQueue<u4> cq(1 MiB);

    // push all at once
    cq.push_back((const u4*)input.data(), N);
    check(cq.size() == N, "size after push");

    // check access
    auto rand_idx = generate_random<u4>(N, 10);
    parlay::sequence<u4> diff(1024);
    parlay::for_each(diff, [&](size_t i){
        diff[i] = input[rand_idx[i]] - cq[rand_idx[i]];
    });
    check(parlay::reduce(diff)==0, "Select unsorted");

    // pop all at once
    size_t n_popped = cq.pop_front(buffer.data(), N + 123456);
    check(n_popped == N, "number of items popped.");
    check(parlay::equal(input, buffer), "push/pop all at once");
    check(cq.size() == 0, "size after all popped");

    p_fill<u4>(buffer.data(), buffer.size(), 0);

    // push one by one
    size_t remaining = N, offset = 0;
    while (remaining) {
        size_t n = random() % (2 MiB);
        n = MIN(n, remaining);
        cq.push_back(input.data() + offset, n);
        remaining -= n;
        offset += n;
    }
    check(cq.size() == N, "size after push");

    // check access
    rand_idx = generate_random<u4>(N, 1024);
    parlay::for_each(diff, [&](size_t i){
        diff[i] = input[rand_idx[i]] - cq[rand_idx[i]];
    });
    check(parlay::reduce(diff)==0, "Select unsorted");

    // pop one by one
    offset = 0;
    size_t num_popped = 0;
    while (true) {
        size_t n = random() % (2 MiB);
        size_t np = cq.pop_front(buffer.data() + offset, n);
        if (!np) break;
        num_popped += np;
        offset += np;
    }
    check(num_popped == N, "number of items popped");
    check(cq.size() == 0, "size after all popped");
    check(parlay::equal(input, buffer), "push/pop in chunks");

}