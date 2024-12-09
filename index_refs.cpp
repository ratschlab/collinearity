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
#include "cqutils.h"

#define SANITY_CHECKS 1
#define BLOCK_SZ (32 MiB)

using namespace klibpp;

#define encode_dna(x) (((x)>>1)&3)

parlay::sequence<u4> create_kmers(const std::string& sequence, int k, int sigma) {
    size_t n = sequence.size();
    if (k > n) return {};
    auto indices = parlay::iota(n - k + 1);
    return parlay::map(indices, [&](size_t i) {
        u4 kmer = 0;
        for (int j = 0; j < k; ++j) kmer = kmer * sigma + encode_dna(sequence[i+j]);
        return kmer;
    });
}

parlay::sequence<u8> create_addresses(u8 ref_id, u8 num_kmers, parlay::sequence<u8> &addresses) {
//    parlay::sequence<u8> addresses(num_kmers, (((u8)ref_id)<<32));
    addresses.resize(num_kmers);
    parlay::for_each(addresses, [&](size_t i) {
        addresses[i] = addresses[i] | ((u8)i);
    });
    return addresses;
}

struct index_t {
    u4 *keys, *keycounts;
    std::vector<compressed_array_t<u8>> values;
};

void * process_fasta(const char* fasta_filename, int k, int sigma) {
    std::vector<std::string> headers;
    // read sequences and convert to key-value pairs
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    auto ks = make_kstream(fd, read, mode::in);

    CQueue<u4> q_keys(BLOCK_SZ, false);
    CQueue<u8> q_values(BLOCK_SZ, true);
    size_t total_nk = 0, total_nbytes = 0;
    uint64_t ref_id = 0;
    parlay::sequence<u8> addresses;
    while (ks >> record) {
        if (record.seq.size() > k) {
            headers.push_back(record.name);
            auto kmers = create_kmers(record.seq, k, sigma);
            const size_t nk = kmers.size();
            q_keys.push_back(kmers.data(), nk);
            create_addresses(ref_id, nk, addresses);
            q_values.push_back(addresses.data(), nk);
            total_nk += nk;
            ref_id++;
            sitrep("Generated %zd tuples from %zd sequences", total_nk, ref_id);
        }
    }
    close(fd);

    printf("\n");
    info("Generated %zd tuples from %zd sequences", total_nk, ref_id);

    auto buf = malloc(BLOCK_SZ * 12);
    cq_sort_by_key(q_keys, q_values, BLOCK_SZ, buf);
    info("Sorted.");

    // compress by q_keys
    CQueue<u4> q_unique_keys(BLOCK_SZ, true);
    CQueue<u4> q_counts(BLOCK_SZ, false);
    cq_count_unique(q_keys, BLOCK_SZ, buf, q_unique_keys, q_counts);
    free(buf);
    info("Counted unique keys.");

    // create array of q_values for each key
    auto index = new index_t;
    size_t nk = q_unique_keys.size();
    expect(nk == q_counts.size());
    index->keys = new u4[nk];
    index->keycounts = new u4[nk];
    q_unique_keys.pop_front(index->keys, nk);
    q_counts.pop_front(index->keycounts, nk);
    size_t max_key_count = *parlay::max_element(parlay::slice(index->keycounts, index->keycounts + nk));
    std::vector<u8> val_buf(max_key_count);
    for (int i = 0; i < nk; ++i) {
        size_t nv = q_values.pop_front(val_buf.data(), index->keycounts[i]);
        verify(nv == index->keycounts[i]);
        val_buf.resize(nv);
        verify(parlay::is_sorted(val_buf));
        index->values.emplace_back(val_buf.data(), nv);
    }
    info("Consolidated values for each unique key.");
    return index;
}

