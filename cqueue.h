//
// Created by Sayan Goswami on 28.11.2024.
//

#ifndef COLLINEARITY_CQUEUE_H
#define COLLINEARITY_CQUEUE_H

#include <deque>
#include "prelude.h"
#include "compressed_array.h"

template <typename T>
struct cqueue_t {
    struct block_t {
        std::vector<T> data;
        size_t start = 0, end = 0;
        block_t(const block_t&) = delete;
        block_t& operator = (const block_t&) = delete;
        block_t(block_t &&other) = default;
        explicit block_t(size_t B) { data.resize(B); }

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
    const size_t B;
    inline size_t size() { return _size; }
    explicit cqueue_t(size_t block_sz): B(block_sz) {}
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
            if (blocks.empty() || !blocks.back().is_pushable()) blocks.emplace_back(B);
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
        if (i < _size) return blocks[i / B].data[i % B];
        else error("Array index out of bounds.");
    }

    T operator[](const size_t i) {
        if (i < _size) return blocks[i / B].data[i % B];
        else error("Array index out of bounds.");
    }
};

template <typename T>
struct CQueue {
    struct buffer_t {
        std::vector<T> data;
        const size_t alloc;
        size_t start, end;

        buffer_t(const buffer_t&) = delete;
        buffer_t& operator = (const buffer_t&) = delete;
        explicit buffer_t(size_t size): data(std::move(std::vector<T>(size))), start(0), end(0), alloc(size) {}
        buffer_t(buffer_t &&other) noexcept : data(std::move(other.data)), start(other.start), end(other.end), alloc(other.alloc) {}

        inline size_t write_to(T* src, size_t n_elements) {
            size_t n = MIN(writeable_size(), n_elements);
            memcpy(data.data() + end, src, n * sizeof(T));
            end += n;
            return n;
        }
        inline size_t read_from(T* dst, size_t n_elements) {
            size_t n = MIN(end - start, n_elements);
            memcpy(dst, data.data() + start, n * sizeof(T));
            start += n;
            if (!can_read_from()) reset();
            return n;
        }
        inline size_t writeable_size() { return alloc - end; }
        inline size_t readable_size() { return end - start; }
        inline void reset() { start = end = 0; }
        inline bool can_write_to() { return writeable_size() > 0; }
        inline bool can_read_from() { return readable_size() > 0; }
    };
private:
    size_t _size;
    buffer_t rb, wb;
    bool sorted;
    std::deque<compressed_array_t<T>> blocks;

public:
    inline void set_sorted(bool is_sorted) { sorted = is_sorted; }
    inline size_t size() { return _size; }
    const size_t CQ_MAXCOUNT;
    CQueue(const CQueue<T>&) = delete;
    CQueue<T>& operator = (const CQueue<T>&) = delete;
    CQueue(CQueue<T> &&other) noexcept :
        CQ_MAXCOUNT(other.CQ_MAXCOUNT), _size(other._size), sorted(other.sorted), rb(std::move(other.rb)), wb(std::move(other.wb))  {
        blocks = std::move(other.blocks);
    }
    explicit CQueue(size_t blocksize = (32 MiB), bool sorted=false) :
        CQ_MAXCOUNT(blocksize), _size(0), rb(buffer_t(blocksize)), wb(buffer_t(blocksize)), sorted(sorted) {}


    void clear() {
        blocks.clear();
        rb.reset(), wb.reset();
        _size = 0;
    }

    /**
     * Bulk-push elements into the queue.
     * @param src - source
     * @param n_elements - number of elements
     */
    void push_back(const T* src, const size_t n_elements) {
        auto arr = (T*)src;
        size_t n_remaining = n_elements;
        while (n_remaining) {
            if (!wb.can_write_to()) {
                blocks.emplace_back(wb.data.data(), wb.end, sorted);
                expect(blocks.back().n == CQ_MAXCOUNT);
                wb.reset();
            }
            size_t n = wb.write_to(arr, n_remaining);
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
            if (rb.can_read_from()) {   // I already have some data in the read buffer
                size_t n = rb.read_from(dst, remaining);
                dst += n, remaining -= n, _size -= n;
            } else if (!blocks.empty()) {   // I have consumed the read buffer, decompress some more
                auto &block = blocks.front();
                block.decompress(rb.data.data());
                rb.end = CQ_MAXCOUNT;
                blocks.pop_front();
            } else if (wb.can_read_from()) {    // consumed all compressed blocks. process my last uncompressed write block if there
                size_t n = wb.read_from(dst, remaining);
                dst += n, remaining -= n, _size -= n;
            } else break;   // I have consumed all data
        }
        return n_elements - remaining;
    }

    /**
     * Pop n elements from this queue and push into another queue
     * @param dst destination queue
     * @param n_elements number of elements be popped
     * @return number of elements actually popped
     */
    size_t pop_front(CQueue<T> &dst, size_t n_elements) {
        auto buffer = new T[CQ_MAXCOUNT];
        size_t remaining = n_elements;
        while (remaining) {
            size_t n = MIN(remaining, CQ_MAXCOUNT);
            size_t n1 = pop_front(buffer, n);
            if (!n1) break;
            dst.push_back(buffer, n1);
            remaining -= n1;
        }
        return n_elements - remaining;
    }

    /**
     * Get the i-th element
     * @param i index
     * @return the element at index i
     */
    const T& operator[](const size_t i) const {
        if (i < _size) {
            if (i < (blocks.size() * CQ_MAXCOUNT)) return blocks.at(i / CQ_MAXCOUNT)[(i % CQ_MAXCOUNT)];
            else return wb.data[i % CQ_MAXCOUNT];
        }
        else error("Array index out of bounds.");
    }

    T operator[](const size_t i) {
        if (i < _size) {
            if (i < (blocks.size() * CQ_MAXCOUNT)) return blocks.at(i / CQ_MAXCOUNT)[(i % CQ_MAXCOUNT)];
            else return wb.data[i % CQ_MAXCOUNT];
        }
        else error("Array index out of bounds.");
    }

//    void copyTo(CQueue<T, sorted> *dst) {
//        for (auto &block : blocks) {
//            dst->push_back(block.data, block.n);
//        }
//    }
//
//    void copyTo(T *dst) {
//        size_t off = 0;
//        for (auto &block : blocks) {
//            memcpy(dst + off, block.data, block.n * sizeof(T));
//            off += block.n;
//        }
//    }
};

#endif //COLLINEARITY_CQUEUE_H
