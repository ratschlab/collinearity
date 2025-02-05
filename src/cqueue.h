//
// Created by Sayan Goswami on 28.11.2024.
//

#ifndef COLLINEARITY_CQUEUE_H
#define COLLINEARITY_CQUEUE_H

#include <deque>
#include "prelude.h"
#include "compressed_array.h"

#define BLOCK_SZ (256 MiB)

template <typename T>
struct cqueue_t {
    struct block_t {
        std::vector<T> data;
        size_t start = 0, end = 0;
        block_t(const block_t&) = delete;
        block_t& operator = (const block_t&) = delete;
        block_t(block_t &&other) = default;
        block_t() { data.resize(BLOCK_SZ); }

        /** how many items can I push to the buffer */
        inline size_t pushable_size() { return data.size() - end; }
        /** how many items can I pop from the buffer */
        inline size_t poppable_size() { return end - start; }
        /** can I push items into the block? */
        inline bool is_pushable() { return pushable_size(); }
        /** can I pop items from the block? */
        inline bool is_poppable() { return poppable_size(); }
        /** reset buffer start and end */
        inline void reset() { start = end = 0; }
        /**
         * Push up to `n_elements` items into the buffer
         * @param src where to push from
         * @param n_elements max. number of items I'll try to push
         * @return the number of items I could actually push based on the available space
         */
        inline size_t push(T* src, size_t n_elements) {
            size_t n = MIN(pushable_size(), n_elements);
            memcpy(data.data() + end, src, n * sizeof(T));
            end += n;
            return n;
        }
        /**
         * Pop up to `n_elements` items from the buffer
         * @param dst where to pop into
         * @param n_elements max. number of items I'll try to pop
         * @return the number of items I could actually pop based on the number of items I have
         */
        inline size_t pop(T* dst, size_t n_elements) {
            size_t n = MIN(end - start, n_elements);
            memcpy(dst, data.data() + start, n * sizeof(T));
            start += n;
            if (!poppable_size()) reset();
            return n;
        }
    };
private:
    size_t _size = 0;
    std::deque<block_t> blocks;
public:
    inline size_t size() { return _size; }
    cqueue_t() = default;
    cqueue_t(const cqueue_t<T>&) = delete;
    cqueue_t<T>& operator = (const cqueue_t<T>&) = delete;
    cqueue_t(cqueue_t<T> &&other) = default;

    /**
     * Bulk-push elements into the queue.
     * @param src - source
     * @param n_elements - number of elements
     */
    void push_back(const T* src, const size_t n_elements) {
        auto arr = (T*)src;
        size_t n_remaining = n_elements;
        while (n_remaining) {
            if (blocks.empty() || !blocks.back().is_pushable()) blocks.emplace_back();
            size_t n = blocks.back().push(arr, n_remaining);
            arr += n, n_remaining -= n, _size += n;
        }
    }

    /**
     * Bulk-pop elements from the queue.
     * @param dst destination array
     * @param n_elements number of elements requested be popped
     * @return the number of items actually popped
     */
    size_t pop_front(T* dst, const size_t n_elements) {
        size_t remaining = n_elements;
        while (remaining) {
            if (!blocks.empty() && !blocks.front().is_poppable()) blocks.pop_front();
            if (blocks.empty()) break;
            auto &block = blocks.front();
            if (block.is_poppable()) {
                size_t n = block.pop(dst, remaining);
                dst += n, remaining -= n, _size -= n;
            }
        }
        return n_elements - remaining;
    }

    /**
     * Pop n elements from this queue and push into another queue
     * @param dst destination queue
     * @param n_elements number of elements be popped
     * @return number of elements actually popped
     */
    size_t pop_front(cqueue_t &dst, size_t n_elements) {
        size_t remaining = n_elements;
        while (remaining) {
            if (!blocks.empty() && !blocks.front().is_poppable()) blocks.pop_front();
            if (blocks.empty()) break;
            auto &block = blocks.front();
            if (block.is_poppable()) {
                size_t n = MIN(block.poppable_size(), remaining);
                dst.push_back(block.data.data() + block.start, n);
                block.start += n, remaining -= n, _size -= n;
            }
        }
        return n_elements - remaining;
    }

    /**
     * clear the blocks
     */
    void clear() {
        blocks.clear();
        _size = 0;
    }

    /**
     * Get the i-th element
     * @param i index
     * @return the element at index i
     */
    const T& operator[](const size_t i) const {
        if (i < _size) return blocks[i / BLOCK_SZ].data[i % BLOCK_SZ];
        else log_error("Array index out of bounds.");
    }

    T operator[](const size_t i) {
        if (i < _size) return blocks[i / BLOCK_SZ].data[i % BLOCK_SZ];
        else log_error("Array index out of bounds.");
    }
};

#endif //COLLINEARITY_CQUEUE_H
