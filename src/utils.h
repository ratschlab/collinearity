//
// Created by Sayan Goswami on 28.01.2025.
//

#ifndef COLLINEARITY_UTILS_H
#define COLLINEARITY_UTILS_H

#include "prelude.h"

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
        return encode_kmer(sequence.data() + i, k, sigma, encoder);
    });
}

static inline parlay::sequence<u8> create_addresses(u8 id, u8 num_kmers) {
    return parlay::tabulate(num_kmers, [&](size_t i) {
        return make_key_from(id, i);
    });
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
