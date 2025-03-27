//
// Created by Sayan Goswami on 29.11.2024.
//

#ifndef COLLINEARITY_TESTS_H
#define COLLINEARITY_TESTS_H

#include "../src/prelude.h"
#include "../src/parlay_utils.h"

#undef verify

#define pass(fmt, ...) fprintf(stderr, "" GRN "[PASS]: " fmt RESET "\n", ##__VA_ARGS__)
#define fail(fmt, ...) fprintf(stderr, "" RED "[FAIL]: " fmt RESET "\n", ##__VA_ARGS__)
#define check(expression, description) ((expression)? pass(description) : fail(description))
#define verify(expression) (expression)? pass("" #expression "") : fail(""  #expression  "")

#define PRINT_VECTOR(v, fmt) ({ \
    for (const auto&x : v) printf(fmt, x); \
    printf("\n"); \
})

#define PRINT_VECTOR1(v, fmt, ...) ({ \
    for (const auto& x : v) printf(fmt, ##__VA_ARGS__); \
    printf("\n"); \
})

template <typename T>
static parlay::sequence<T> generate_random(size_t max, size_t count, int seed=0) {
    parlay::random_generator gen(seed);
    std::uniform_int_distribution<T> dis(0, max);
    auto result = parlay::tabulate(count, [&](size_t i) {
        auto r = gen[i];
        return dis(r);
    });
    return result;
}

#define f(x) void f##x()

f(a0);
f(a1);
f(a2);
f(a3);

#endif //COLLINEARITY_TESTS_H
