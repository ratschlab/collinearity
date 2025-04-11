//
// Created by Sayan Goswami on 28.01.2025.
//

#ifndef COLLINEARITY_UTILS_H
#define COLLINEARITY_UTILS_H

#include "prelude.h"

#define LOW32(x) ((x) & 0xffffffff)
#define HIGH32(x) (((x)>>32) & 0xffffffff)
#define MAKE64(hi, lo) ((((u8)(hi))<<32) | (lo))

static auto encode_dna = [](char x) { return (u4)(((x) >> 1) & 3); };
static auto encode_qsig = [](u1 x) { return x; };

static u8 ipow(u8 base, u8 exp) {
    u8 result = 1;
    for (;;) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }
    return result;
}

template<typename T>
static inline  size_t upper_bound(const T *arr, size_t start, size_t end, const T& key) {
    size_t i;
    while (start < end) {
        i = (start + end) / 2;
        if (key < arr[i]) end = i;
        else start = i + 1;
    }
    return start;
}

template<typename T>
static inline  size_t lower_bound(const T *arr, size_t start, size_t end, const T& key) {
    size_t i;
    while (start < end) {
        i = (start + end) / 2;
        if (key > arr[i]) start = i + 1;
        else end = i;
    }
    return start;
}

template <typename T, typename Encoder>
static inline u4 encode_kmer(const T& seq, int k, int sigma, Encoder encoder) {
    u8 v = 0;
    for (u4 i = 0; i < k; ++i) v = v * sigma + encoder(seq[i]);
    return v;
}

template <typename T, typename Encoder>
static inline parlay::sequence<u4> create_kmers(const T& sequence, int k, int sigma, Encoder encoder) {
    size_t n = sequence.size();
    if (k > n) return {};
    return parlay::tabulate(n - k + 1, [&](size_t i){
        return encode_kmer(sequence.begin() + i, k, sigma, encoder);
    });
}

template <typename T, typename Encoder>
static inline parlay::sequence<u4> create_kmers_1t(const T& sequence, int k, int sigma, Encoder encoder) {
    const u8 M = ipow(sigma, k-1);
    const u4 n = sequence.size(), n_keys = n - k + 1;
    expect(n > k);
    parlay::sequence<u4> keys(n_keys);
    keys[0] = encode_kmer(sequence, k, sigma, encoder);
    for (u4 i = k, j = 1; i < n; ++i, ++j) {
        keys[j] = (keys[j-1] - encoder(sequence[j-1]) * M) * sigma + encoder(sequence[i]);
    }
    return keys;
}

static inline parlay::sequence<char> revcmp(std::string &seq) {
    const auto n = seq.size();
    parlay::sequence<char> revseq(n);
    for (int i = 0; i < n; ++i) revseq[i] = "TGAC"[(seq[n-1-i] >> 1) & 3];
    return revseq;
}

static inline parlay::sequence<char> revcmp_par(std::string &seq) {
    const auto n = seq.size();
    return parlay::tabulate(n, [&](size_t i){
        return "TGAC"[(seq[n-1-i] >> 1) & 3];
    });
}

template <template<typename> class Sequence, typename T>
static void prettyPrintVector(const Sequence<T>& vec, const char *fmt, size_t threshold = 20) {
    size_t n = vec.size();
    printf("[");
    if (n <= threshold) {
        for (size_t i = 0; i < n-1; ++i) {
            printf(fmt, vec[i]);
            printf(", ");
        }
        printf(fmt, vec[n-1]);
    } else {
        size_t headCount = threshold / 2; // Elements to show from the head
        size_t tailCount = threshold / 2;
        for (size_t i = 0; i < headCount; ++i) {
            printf(fmt, vec[i]);
            printf(", ");
        }
        printf("... (%d more elements) ...", n - threshold);
        for (size_t i = n - tailCount; i < n-1; ++i) {
            printf(fmt, vec[i]);
            printf(", ");
        }
        printf(fmt, vec[n-1]);
    }
    printf("]\n");
    fflush(stdout);
}

// Pretty-print a vector, showing head and tail for long vectors
template <template<typename> class Sequence, typename T>
static void prettyPrintVector(const Sequence<T>& vec, size_t threshold = 20) {
    int precision = 5;
    size_t n = vec.size();

    std::cout << "[";
    if (n <= threshold) {
        // Print the full vector if it's small enough
        for (size_t i = 0; i < n; ++i) {
            if constexpr (std::is_floating_point_v<T>) {
                std::cout << std::fixed << std::setprecision(precision);
            }
            std::cout << vec[i];
            if (i != n - 1) {
                std::cout << ", ";
            }
        }
    } else {
        // Print head, ellipsis, and tail for long vectors
        size_t headCount = threshold / 2; // Elements to show from the head
        size_t tailCount = threshold / 2; // Elements to show from the tail

        for (size_t i = 0; i < headCount; ++i) {
            if constexpr (std::is_floating_point_v<T>) {
                std::cout << std::fixed << std::setprecision(precision);
            }
            std::cout << vec[i] << ", ";
        }

        std::cout << "... (" << n - threshold << " more elements) ...";

        for (size_t i = n - tailCount; i < n; ++i) {
            if constexpr (std::is_floating_point_v<T>) {
                std::cout << std::fixed << std::setprecision(precision);
            }
            std::cout << (i == n - tailCount ? " " : ", ") << vec[i];
        }
    }
    std::cout << "]" << std::endl;
}



#endif //COLLINEARITY_UTILS_H
