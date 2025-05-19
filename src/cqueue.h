//
// Created by Sayan Goswami on 28.11.2024.
//

#ifndef COLLINEARITY_CQUEUE_H
#define COLLINEARITY_CQUEUE_H

#include <deque>
#include "prelude.h"
#include "mempool.h"

template <typename T>
struct cqueue_t {
    struct block_t {
        mempool_t<T> &mp;
        T* data = nullptr;
        size_t start = 0, end = 0;
        block_t(const block_t&) = delete;
        block_t& operator = (const block_t&) = delete;
        block_t(block_t &&other) = default;
        block_t() : mp(mempool_t<T>::getInstance()), data(mp.reserve()) { }
        ~block_t() {
            if (data) mp.release(data);
            data = nullptr;
        }

        /** how many items can I push to the buffer */
        inline size_t pushable_size() { return MEMPOOL_BLOCKSZ - end; }
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
            memcpy(data + end, src, n * sizeof(T));
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
            memcpy(dst, data + start, n * sizeof(T));
            start += n;
            if (!poppable_size()) reset();
            return n;
        }
    };
private:
    size_t _size = 0;
    std::deque<block_t> blocks;
public:
    typedef T value_type;
    inline size_t size() const { return _size; }
    inline bool empty() const { return _size == 0; }
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
            if (blocks.empty()) break;
            auto &block = blocks.front();
            if (block.is_poppable()) {
                size_t n = block.pop(dst, remaining);
                dst += n, remaining -= n, _size -= n;
            }
            if (!blocks.empty() && !blocks.front().is_poppable()) blocks.pop_front();
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
            if (blocks.empty()) break;
            auto &block = blocks.front();
            if (block.is_poppable()) {
                size_t n = MIN(block.poppable_size(), remaining);
                dst.push_back(block.data + block.start, n);
                block.start += n, remaining -= n, _size -= n;
            }
            if (!blocks.empty() && !blocks.front().is_poppable()) blocks.pop_front();
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
        if (i < _size) return blocks[MPDIV(i)].data[MPMOD(i)];
        else log_error("Array index out of bounds.");
    }

    T operator[](const size_t i) {
        if (i < _size) return blocks[MPDIV(i)].data[MPMOD(i)];
        else log_error("Array index out of bounds.");
    }

    class const_iterator {
    private:
        const cqueue_t<T> *q;
        size_t i;
    public:
        using value_type = T;
        using difference_type = size_t;
        using pointer = T*;
        using reference = T&;
        using iterator_category = std::forward_iterator_tag;

        // Constructor
        explicit const_iterator(const cqueue_t<unsigned int> *q) : q(q), i(0) {}
        const_iterator(const cqueue_t<T> *q, size_t i) : q(q), i(i) {}

        // Dereference operator
        T operator*() { return (*q)[i]; }

        // Arrow operator
        T* operator->() { return &((*q)[i]); }

        // Prefix increment
        const_iterator& operator++() {
            ++i;
            return *this;
        }

        // Postfix increment
        const_iterator operator++(int) {
            const_iterator temp = *this;
            ++(*this);
            return temp;
        }

        // Difference between two iterators
        difference_type operator-(const const_iterator& other) const {
            return i - other.i;
        }

        // Subscript operator
        T& operator[](difference_type n) const {
            return *((*q)[i + n]);
        }

        // Comparison operators
        bool operator==(const const_iterator& other) const { return i == other.i; }
        bool operator!=(const const_iterator& other) const { return i != other.i; }
        bool operator<(const const_iterator& other) const { return i < other.i; }
        bool operator<=(const const_iterator& other) const { return i <= other.i; }
        bool operator>(const const_iterator& other) const { return i > other.i; }
        bool operator>=(const const_iterator& other) const { return i >= other.i; }
    };

    // Begin and End functions
    const_iterator begin() const { return const_iterator(this); }
    const_iterator end() const { return const_iterator(this, size()); }
    const_iterator it(size_t i) const { return const_iterator(this, i); }

    void dump(std::ofstream &fs) {
        dump_values(fs, _size);
        for (auto &block : blocks)
            dump_data(fs, block.data, block.end - block.start);
    }

    void load(std::ifstream &fs) {
        if (_size) log_error("Loading a file into a non-empty queue is not supported.");
        size_t size;
        load_values(fs, &size);
        _size = size;
        while (size) {
            size_t rc = std::min(MEMPOOL_BLOCKSZ, size);
            blocks.emplace_back();
            load_data(fs, blocks.back().data, rc);
            blocks.back().end = rc;
            size -= rc;
        }
    }
};

#endif //COLLINEARITY_CQUEUE_H
