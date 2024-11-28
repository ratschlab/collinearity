//
// Created by Sayan Goswami on 27.11.2024.
//

#ifndef COLLINEARITY_COMPRESSED_ARRAY_H
#define COLLINEARITY_COMPRESSED_ARRAY_H

#include <cstring>
#include "vbyte.h"
#include "prelude.h"

template <typename T, bool sorted>
struct compressed_array_t {
    void *data;
    const size_t n;
    size_t nc;
    bool compressed;
    compressed_array_t(T *arr, size_t n): n(n) {
        if constexpr (std::is_same<T, u4>::value) {
            if constexpr (sorted) {
                nc = vbyte_compressed_size_sorted32(arr, n, 0);
            } else {
                nc = vbyte_compressed_size_unsorted32(arr, n);
            }
        } else if constexpr (std::is_same<T, u8>::value) {
            if constexpr (sorted) {
                nc = vbyte_compressed_size_sorted64(arr, n, 0);
            } else {
                nc = vbyte_compressed_size_unsorted64(arr, n);
            }
        } else {
            error("Type has to be either u4 or u8.");
        }

        if (nc < n * sizeof(T)) {
            data = malloc(nc);
            size_t nc1;

            if constexpr (std::is_same<T, u4>::value) {
                if constexpr (sorted) {
                    nc1 = vbyte_compress_sorted32(arr, (u1*)data, n, 0);
                } else {
                    nc1 = vbyte_compress_unsorted32(arr, (u1*)data, n);
                }
            } else if constexpr (std::is_same<T, u8>::value) {
                if constexpr (sorted) {
                    nc1 = vbyte_compress_sorted64(arr, (u1*)data, n, 0);
                } else {
                    nc1 = vbyte_compress_unsorted64(arr, (u1*)data, n);
                }
            } else {
                error("Type has to be either u4 or u8.");
            }
            if(nc1 != nc)
                warn("Calculated size is not equal to actual size.");
            compressed = true;
        } else {
            nc = n * sizeof(T);
            data = malloc(nc);
            memcpy(data, arr, nc);
            compressed = false;
        }
    }

    void decompress(T* out, const size_t size) {
        expect(size >= n);
        if (!compressed) memcpy(out, data, nc);
        else {
            if constexpr (std::is_same<T, u4>::value) {
                if constexpr (sorted) {
                    vbyte_uncompress_sorted32((u1*)data, out, 0, n);
                } else {
                    vbyte_uncompress_unsorted32((u1*)data, out, n);
                }
            } else if constexpr (std::is_same<T, u8>::value) {
                if constexpr (sorted) {
                    vbyte_uncompress_sorted64((u1*)data, out, 0, n);
                } else {
                    vbyte_uncompress_unsorted64((u1*)data, out, n);
                }
            } else {
                error("Type has to be either u4 or u8.");
            }
        }
    }

    ~compressed_array_t() {
        free(data);
    }
};

#endif //COLLINEARITY_COMPRESSED_ARRAY_H
