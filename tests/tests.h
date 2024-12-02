//
// Created by Sayan Goswami on 29.11.2024.
//

#ifndef COLLINEARITY_TESTS_H
#define COLLINEARITY_TESTS_H

#include "../prelude.h"
#include "../parlay_utils.h"

#define pass(fmt, ...) fprintf(stderr, "" GRN "[PASS]: " fmt RESET "\n", ##__VA_ARGS__)
#define fail(fmt, ...) fprintf(stderr, "" RED "[FAIL]: " fmt RESET "\n", ##__VA_ARGS__)
#define check(expression, description) ((expression)? pass(description) : fail(description))

template <typename T>
static parlay::sequence<T> generate_random(size_t max, size_t count) {
    parlay::random_generator gen;
    std::uniform_int_distribution<T> dis(0, max);
    return parlay::tabulate(count, [&](size_t i) {
        auto r = gen[i];
        return dis(r);
    });
}

void test_compressed_array();
void test_cqueue();
void test_cqutils();

#endif //COLLINEARITY_TESTS_H
