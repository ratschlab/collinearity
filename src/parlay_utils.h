//
// Created by Sayan Goswami on 28.11.2024.
//

#ifndef COLLINEARITY_PARLAY_UTILS_H
#define COLLINEARITY_PARLAY_UTILS_H

#include "prelude.h"
#include "parlay/primitives.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/io.h"

/**
 * Fill an array with value
 * @tparam T type
 * @param arr array
 * @param count size of array
 * @param value value to fill with
 */
template <typename T>
static inline void p_fill(T *arr, size_t count, const T& value) {
    parlay::parallel_for(0, count, [&](size_t i){
        arr[i] = value;
    });
}

/**
 * Create a sequence and store it into an array
 * @tparam T type
 * @param arr array
 * @param count length of array
 * @param start start of sequence
 * @param step step of sequence
 */
template <typename T>
static inline void p_seq(T *arr, size_t count, const T& start=0, const T& step=1) {
    parlay::parallel_for(0, count, [&](size_t i){
        arr[i] = start + i*step;
    });
}

/**
 * Scatter elements from `in` to `out` according to `index` such that `out[index[i]] = in[i]`
 * @tparam T array type
 * @tparam I index type
 * @param in input array
 * @param n size of `in`
 * @param index index used to map elements from in to out
 * @param out the output array. can be same as in but otherwise should not overlap with it
 */
template <typename T, typename I>
static inline void p_scatter(T *in, size_t n, I *index, T *out) {
    parlay::sequence<T> tmp;
    bool in_place = (in == out);
    if (in_place) {
        tmp = parlay::sequence<T>(n);
        out = tmp.data();
    }
    parlay::parallel_for(0, n, [&](size_t i){ out[index[i]] = in[i]; });
    if (in_place) memcpy(in, out, n * sizeof(T));
}

/**
 * Sort tuples by one of its elements
 * @tparam T1 key type
 * @tparam T2 value type
 * @param a key array
 * @param b value array
 * @param n length of the input arrays
 */
template <typename T1, typename T2>
static inline void sort_by_key(T1* a, T2 *b, size_t n) {
    auto index = parlay::rank(parlay::slice(a, a + n));
    p_scatter(a, n, index.data(), a);
    p_scatter(b, n, index.data(), b);
}

template <typename T>
static inline bool is_sorted(T *a, size_t n) {
    return parlay::is_sorted(parlay::slice(a, a + n));
}

/**
 * Sort tuples by one of the elements
 * @tparam T1 key type
 * @tparam T2 value type 1
 * @tparam T3 value type 2
 * @param a key array
 * @param b value array 1
 * @param c value array 2
 * @param n size of the unput arrays
 */
template <typename T1, typename T2, typename T3>
static inline void sort_by_key(T1* a, T2 *b, T3 *c, size_t n) {
    auto index = parlay::rank(parlay::slice(a, a + n));
    p_scatter(a, n, index.data(), a);
    p_scatter(b, n, index.data(), b);
    p_scatter(c, n, index.data(), c);
}

#endif //COLLINEARITY_PARLAY_UTILS_H
