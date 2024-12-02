//
// Created by Sayan Goswami on 27.11.2024.
//

#include <vector>
#include "prelude.h"
#include "kseq++/kseq++.hpp"
#include "vbyte.h"
#include <type_traits>
#include "compressed_array.h"
#include "parlay_utils.h"

#define SANITY_CHECKS 1
#define BLOCK_SZ (32 MiB)

using namespace klibpp;

#define encode_dna(x) (((x)>>1)&3)

parlay::sequence<u4> create_kmers(const std::string& sequence, int k, int sigma) {
    size_t n = sequence.size();
    if (k > n) return {};
    auto indices = parlay::iota(n - k + 1);
    auto kmers = parlay::map(indices, [&](size_t i) {
        u4 kmer = 0;
        for (int j = 0; j < k; ++j) kmer = kmer * sigma + encode_dna(sequence[i+j]);
        return kmer;
    });
    return kmers;
}

struct kmer_chunk_t {
    const std::string name;
    compressed_array_t<u4> kmers;
    kmer_chunk_t(std::string &name, u4 *data, size_t length): name(name), kmers(data, length) {}
    [[nodiscard]] size_t size() const { return kmers.n; }
    [[nodiscard]] size_t compressed_size() const {return kmers.nc; }
};

struct tuple_chunk_t {
    compressed_array_t<u4> keys;
    compressed_array_t<u4> rids;
    compressed_array_t<u4> poss;
    parlay::sequence<u4> tmp_keys, tmp_rids, tmp_poss;
    const size_t n;

    tuple_chunk_t(u4 *keys, u4 *rids, u4 *poss, size_t n): n(n), keys(keys, n), rids(rids, n), poss(poss, n) {}

    inline void summarize() const {
        size_t csz = keys.nc + rids.nc + poss.nc;
        info("Stored %zd tuples in %.2f GB (compression ratio %.2f%%)", n, csz*1.0/(1 GiB), csz * 100.0 / (n *
                sizeof(u4)));
    }

    void decompress(u4 *keyout, u4 *ridout, u4 *posout) {
        keys.decompress(keyout);
        rids.decompress(ridout);
        poss.decompress(posout);
    }

    void decompress(parlay::sequence<u4> &keyout, parlay::sequence<u4> &ridout, parlay::sequence<u4> &posout) {
        keyout.resize(n);
        ridout.resize(n);
        posout.resize(n);
        decompress(keyout.data(), ridout.data(), posout.data());
    }
};

struct sorted_tuples_t {
    std::list<tuple_chunk_t> tuple_chunks;
    parlay::sequence<u4> keys, rids, poss;
    u4 last_key;
    size_t start = 0, end = 0;

    explicit sorted_tuples_t(tuple_chunk_t &&tuple_chunk) { tuple_chunks.push_back(tuple_chunk); }
    explicit sorted_tuples_t(std::list<tuple_chunk_t> &tuple_chunks) {
        while (!tuple_chunks.empty()) {
            auto &&tuples = tuple_chunks.front();
            this->tuple_chunks.push_back(tuples);
            tuple_chunks.pop_front();
        }
    }

    size_t stage() {
        if (start == end) {
            if (tuple_chunks.empty()) return 0;
            auto &tuples = tuple_chunks.front();
            tuples.decompress(keys, rids, poss);
            end = keys.size();
            tuple_chunks.pop_front();
            last_key = keys[end-1];
        }
        return end - start;
    }
    inline size_t consume(size_t n) {
        n = MIN(n, (end-start));
        start += n;
        return n;
    }
    inline void rewind(size_t n) { start -= n; }
};

static void merge_two_sorted_tuples(sorted_tuples_t &a, sorted_tuples_t &b) {
    size_t na = a.stage(), nb = b.stage();
    while (na && nb) {
        // find boundary
        u4 min_key = MIN(a.last_key, b.last_key);
        size_t endA = lower_bound(a.keys.data(), a.start, a.end, min_key);
        size_t endB = lower_bound(b.keys.data(), b.start, b.end, min_key);

    }
}

static void merge_many_sorted_tuples(std::list<tuple_chunk_t> &sorted_tuple_chunks) {
    std::list<sorted_tuples_t> sorted_lists;
    while (!sorted_tuple_chunks.empty()) {
        sorted_lists.emplace_back(sorted_tuple_chunks);
        sorted_tuple_chunks.pop_front();
    }
    while (sorted_lists.size() > 1) {
        auto &&A = sorted_lists.front();
        auto &&B = sorted_lists.front();
        sorted_lists.pop_front();
        sorted_lists.pop_front();

    }
}

static void sort_tuples(parlay::sequence<u4> &keys, parlay::sequence<u4> &rids, parlay::sequence<u4> &poss, size_t n,
                        std::list<tuple_chunk_t> &sorted_tuple_chunks) {
    keys.resize(n);
    rids.resize(n);
    poss.resize(n);
    info("Sorting by position..");
    sort_by_key(poss.data(), keys.data(), rids.data(), n);
    info("Sorting by seq id..");
    sort_by_key(rids.data(), keys.data(), poss.data(), n);
    info("Sorting by keys..");
    sort_by_key(keys.data(), rids.data(), poss.data(), n);

    if (SANITY_CHECKS) {
        info("Checking sanity..");
        verify(parlay::is_sorted(keys));
        auto histogram_keys = parlay::histogram_by_key(keys);
        size_t key_offset = 0, rid_offset = 0;
        for (auto &[key, count]: histogram_keys) {
            auto seq = parlay::sequence<u4>(rids.begin() + key_offset, rids.begin() + key_offset + count);
            verify(parlay::is_sorted(seq));
            auto histogram_rids = parlay::histogram_by_key(seq);
            for (auto &[rid, count] : histogram_rids) {
                auto seq = parlay::sequence<u4>(poss.begin() + rid_offset, poss.begin() + rid_offset + count);
                verify(parlay::is_sorted(seq));
                rid_offset += count;
            }
            key_offset += count;
        }
        info("All good.");
    }

    sorted_tuple_chunks.emplace_back(keys.data(), rids.data(), poss.data(), n);
    sorted_tuple_chunks.back().summarize();
}

std::list<tuple_chunk_t> sort_tuple_blocks(std::list<kmer_chunk_t> &kmer_chunks) {
    std::list<tuple_chunk_t> sorted_tuple_chunks;
    size_t max_buf_sz = BLOCK_SZ;
    parlay::sequence<u4> keys(max_buf_sz);
    parlay::sequence<u4> rids(max_buf_sz);
    parlay::sequence<u4> poss(max_buf_sz);
    size_t remaining = max_buf_sz, offset = 0;
    u4 ref_id = 0;
    while (!kmer_chunks.empty()) {
        auto &kmer_chunk = kmer_chunks.front();
        auto &kmers = kmer_chunk.kmers;
        if (kmers.n > remaining) {
            if (offset == 0) {
                max_buf_sz = kmers.n;
                keys.resize(max_buf_sz);
                kmers.decompress(keys.data());
                p_fill(rids.data(), kmers.n, ref_id);
                p_seq(poss.data(), kmers.n);
                offset = kmers.n;
                remaining = 0;
            }
            sort_tuples(keys, rids, poss, offset, sorted_tuple_chunks);
            offset = 0;
            remaining = max_buf_sz;
        }
        kmers.decompress(keys.data() + offset);
        remaining -= kmers.n;
        offset += kmers.n;
        kmer_chunks.pop_front();
        ref_id++;
    }
    if (offset) sort_tuples(keys, rids, poss, offset, sorted_tuple_chunks);
    return sorted_tuple_chunks;
}

void process_fasta(const char* fasta_filename, int k, int sigma) {
    // read sequences and convert to key-value pairs
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    auto ks = make_kstream(fd, read, mode::in);

    std::list<kmer_chunk_t> kmer_chunks;
    size_t total_nk = 0, total_nbytes = 0;
    while (ks >> record) {
        auto kmers = create_kmers(record.seq, k, sigma);
        const size_t nk = kmers.size();
        kmer_chunks.emplace_back(record.name, kmers.data(), nk);
        auto nc = kmer_chunks.back().compressed_size();
        sitrep("Generated %zd kmers and stored in %zd bytes (%.2f%%).", nk, nc, nc*100.0/(nk*sizeof(u4)));
        total_nbytes += nc;
        total_nk += nk;
    }
    close(fd);

    printf("\n");
    info("In total, stored %zd kmers in %.2f GB (compression ratio %.2f%%)",
         total_nk, total_nbytes * 1.0/(1 GiB), total_nbytes * 100.0/(total_nk * sizeof(u4)));

    auto sorted_blocks = sort_tuple_blocks(kmer_chunks);
    info("Created %zd sorted blocks.", sorted_blocks.size());

}

