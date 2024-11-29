//
// Created by Sayan Goswami on 29.11.2024.
//

#ifndef COLLINEARITY_TESTS_H
#define COLLINEARITY_TESTS_H

#include "../prelude.h"

#define pass(fmt, ...) fprintf(stderr, "" GRN "[PASS]: " fmt RESET "\n", ##__VA_ARGS__)
#define fail(fmt, ...) fprintf(stderr, "" RED "[FAIL]: " fmt RESET "\n", ##__VA_ARGS__)
#define test(expression, description) ((expression)? pass(description) : fail(description))

void test_compressed_array();

#endif //COLLINEARITY_TESTS_H
