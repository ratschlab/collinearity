//
// Created by Sayan Goswami on 27.11.2024.
//

#include "collinearity.h"

#define SANITY_CHECKS 1
#define BLOCK_SZ (32 MiB)

using namespace klibpp;

void create_addresses(u8 ref_id, u8 num_kmers, parlay::sequence<u8> &addresses) {
//    parlay::sequence<u8> addresses(num_kmers, );
    addresses.resize(num_kmers);
    parlay::for_each(parlay::iota(num_kmers), [&](size_t i) {
        addresses[i] = make_key_from(ref_id, i);
    });
}

std::pair<index_t*, std::vector<std::string>> process_fasta(const char* fasta_filename, int k, int sigma) {
    std::vector<std::string> headers;
    // read sequences and convert to key-value pairs
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    if (fd < 0) error("Could not open %s because %s.", fasta_filename, strerror(errno));
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

    auto index = new index_t(q_unique_keys, q_counts, q_values);
    return std::make_pair(index, std::move(headers));
}

