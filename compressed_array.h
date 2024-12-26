//
// Created by Sayan Goswami on 27.11.2024.
//

#ifndef COLLINEARITY_COMPRESSED_ARRAY_H
#define COLLINEARITY_COMPRESSED_ARRAY_H

#include <cstring>
#include "vbyte.h"
#include "prelude.h"

template <typename T>
struct compressed_array_t {
    size_t n;
    std::vector<T> v;
    compressed_array_t() : n(0) {}
    compressed_array_t(const T *arr, const size_t n, bool sorted = false) : n(n) {
        compress(arr, n, sorted);
    }
    const T& operator[](size_t i) const { return v[i]; }
    T operator[](size_t i) { return v[i]; }

    void compress(const T *arr, const size_t size, bool is_sorted = false) {
        n = size;
        v.assign(arr, arr + n);
    }

    size_t decompress(T* out) {
        std::memcpy(out, v.data(), n * sizeof(T));
        return n;
    }
};

//template <typename T>
//struct compressed_array_t {
//    size_t n;
//    size_t nc;
//    std::vector<u1> v;
//    const bool sorted;
//    bool compressed;
//    compressed_array_t(): n(0), nc(0), sorted(false), compressed(false) {}
//    compressed_array_t(const T *arr, const size_t n, bool sorted = false): n(n), sorted(sorted), nc(0), compressed(false) {
//        compress(arr, n, sorted);
//    }
//
//    void compress(const T *arr, const size_t size, bool is_sorted = false) {
//        n = size;
//        if constexpr (std::is_same<T, u4>::value) {
//            if (is_sorted) {
//                nc = vbyte_compressed_size_sorted32(arr, size, 0);
//            } else {
//                nc = vbyte_compressed_size_unsorted32(arr, size);
//            }
//        } else if constexpr (std::is_same<T, u8>::value) {
//            if (is_sorted) {
//                nc = vbyte_compressed_size_sorted64(arr, size, 0);
//            } else {
//                nc = vbyte_compressed_size_unsorted64(arr, size);
//            }
//        } else {
//            error("Type has to be either u4 or u8.");
//        }
//
//        if (nc < size * sizeof(T)) {
////            data = malloc(nc);
//            v.resize(nc);
//            size_t nc1;
//
//            if constexpr (std::is_same<T, u4>::value) {
//                if (is_sorted) {
//                    nc1 = vbyte_compress_sorted32(arr, (u1*)v.data(), 0, size);
//                } else {
//                    nc1 = vbyte_compress_unsorted32(arr, (u1*)v.data(), size);
//                }
//            } else if constexpr (std::is_same<T, u8>::value) {
//                if (is_sorted) {
//                    nc1 = vbyte_compress_sorted64(arr, (u1*)v.data(), 0, size);
//                } else {
//                    nc1 = vbyte_compress_unsorted64(arr, (u1*)v.data(), size);
//                }
//            } else {
//                error("Type has to be either u4 or u8.");
//            }
//            if(nc1 != nc)
//                warn("Calculated size is not equal to actual size.");
//            compressed = true;
//        } else {
////            warn("Coudn't compress. Compression ratio = %.2f%%", nc * 100.0/(size * sizeof(T)));
//            nc = size * sizeof(T);
////            data = malloc(nc);
//            v.resize(nc);
//            memcpy(v.data(), arr, nc);
//            compressed = false;
//        }
//    }
//
//    size_t decompress(T* out) {
//        if (!compressed) memcpy(out, v.data(), nc);
//        else {
//            if constexpr (std::is_same<T, u4>::value) {
//                if (sorted) {
//                    vbyte_uncompress_sorted32((u1*)v.data(), out, 0, n);
//                } else {
//                    vbyte_uncompress_unsorted32((u1*)v.data(), out, n);
//                }
//            } else if constexpr (std::is_same<T, u8>::value) {
//                if (sorted) {
//                    vbyte_uncompress_sorted64((u1*)v.data(), out, 0, n);
//                } else {
//                    vbyte_uncompress_unsorted64((u1*)v.data(), out, n);
//                }
//            } else {
//                error("Type has to be either u4 or u8.");
//            }
//        }
//        return n;
//    }
//
//    const T& operator[](size_t i) const {
//        if (i < n) {
//            if (compressed) {
//                if constexpr (std::is_same<T, u4>::value) {
//                    if (sorted) {
//                        return vbyte_select_sorted32((u1 *)v.data(), nc, 0, i);
//                    } else {
//                        return vbyte_select_unsorted32((u1 *)v.data(), nc, i);
//                    }
//                } else if constexpr (std::is_same<T, u8>::value) {
//                    if (sorted) {
//                        return vbyte_select_sorted64((u1 *)v.data(), nc, 0, i);
//                    } else {
//                        return vbyte_select_unsorted64((u1 *)v.data(), nc, i);
//                    }
//                } else {
//                    error("Type has to be either u4 or u8.");
//                }
//            } else return ((T*)v.data())[i];
//        }
//        else error("Array index out of bounds.");
//    }
//
//    T operator[](size_t i) {
//        if (i < n) {
//            if (compressed) {
//                if constexpr (std::is_same<T, u4>::value) {
//                    if (sorted) {
//                        return vbyte_select_sorted32((u1 *)v.data(), nc, 0, i);
//                    } else {
//                        return vbyte_select_unsorted32((u1 *)v.data(), nc, i);
//                    }
//                } else if constexpr (std::is_same<T, u8>::value) {
//                    if (sorted) {
//                        return vbyte_select_sorted64((u1 *)v.data(), nc, 0, i);
//                    } else {
//                        return vbyte_select_unsorted64((u1 *)v.data(), nc, i);
//                    }
//                } else {
//                    error("Type has to be either u4 or u8.");
//                }
//            } else return ((T*)v.data())[i];
//        }
//        else error("Array index out of bounds.");
//    }
//};

#endif //COLLINEARITY_COMPRESSED_ARRAY_H
