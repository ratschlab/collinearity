//
// Created by Sayan Goswami on 26.12.2024.
//

#ifndef COLLINEARITY_MEMPOOL_H
#define COLLINEARITY_MEMPOOL_H

#include "prelude.h"
#include <queue>

// our block size is 2^26 = 64Mi
#define BLOCKSIZE_BITS 26
#define MEMPOOL_BLOCKSZ (1UL<<BLOCKSIZE_BITS)
#define MPDIV(x) ((x)>>BLOCKSIZE_BITS)
#define MPMOD(x) ((x) & (MEMPOOL_BLOCKSZ-1))

template <typename T>
class mempool_t {
    int n_reserved = 0, n_total_allocated = 0;
    std::queue<T*> chunks;
    mempool_t() = default;
    ~mempool_t() { shrink(); }
public:
    mempool_t(const mempool_t&) = delete;
    mempool_t& operator=(const mempool_t&) = delete;
    static mempool_t& getInstance() {
        static mempool_t<T> instance;
        return instance;
    }

    T* reserve() {
        T *chunk;
        if (chunks.empty()) {
            chunk = new T[MEMPOOL_BLOCKSZ];
            n_total_allocated++;
        } else {
            chunk = chunks.front();
            chunks.pop();
        }
        n_reserved++;
        return chunk;
    }

    void release(T *chunk) {
        chunks.push(chunk);
        --n_reserved;
    }

    void shrink() {
        while (!chunks.empty()) {
            T *chunk = chunks.front();
            chunks.pop();
            delete[] chunk;
            --n_total_allocated;
        }
    }

    void print_usage() {
        log_debug(LOW, "MP(%zd): %d blocks in use, %d blocks allocated in total.",
                  sizeof(T), n_reserved, n_total_allocated);
    }
};

#define PRINT_MEM_USAGE(T) mempool_t<T>::getInstance().print_usage()
#define MEMPOOL_SHRINK(T) mempool_t<T>::getInstance().shrink()

#endif //COLLINEARITY_MEMPOOL_H
