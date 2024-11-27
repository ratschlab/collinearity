//
// Created by Sayan Goswami on 27.11.2024.
//

#include <vector>
#include "prelude.h"
#include <parlay/primitives.h>
#include <parlay/io.h>
#include "kseq++/kseq++.hpp"
#include "vbyte.h"

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

struct kmer_array_t {
    void *data;
    bool compressed;
    size_t nk, nbytes;
    kmer_array_t(void *data, size_t nk, size_t nbytes):
        data(data), nk(nk), nbytes(nbytes), compressed(nk * sizeof(u4) > nbytes) {}
};

struct tuples_t {
    std::vector<std::string> names;
    std::vector<uint64_t> lengths;
    std::vector<kmer_array_t> kmers;
};

void process_fasta(const char* fasta_filename, int k, int sigma) {
    // read sequences and convert to key-value pairs
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    auto ks = make_kstream(fd, read, mode::in);

    tuples_t tuples;
    size_t total_nk = 0, total_nbytes = 0;
    while (ks >> record) {
        tuples.names.emplace_back(record.name);
        auto kmers = create_kmers(record.seq, k, sigma);
        const size_t nk = kmers.size();
        tuples.lengths.push_back(nk);
        size_t nc = vbyte_compressed_size_unsorted32(kmers.data(), nk);
        if (nc < nk * sizeof(u4)) {
            auto data = new uint8_t [nc];
            if(nc != vbyte_compress_unsorted32(kmers.data(), data, nk))
                warn("Calculated size is not equal to actual size.");
            tuples.kmers.emplace_back(data, nk, nc);
        } else {
            nc = nk * sizeof(u4);
            auto data = new u4[nk];
            memcpy(data, kmers.data(), nc);
            tuples.kmers.emplace_back(data, nk, nc);
        }
        sitrep("Generated %zd kmers and stored in %zd bytes (%.2f%%).", nk, nc, nc*100.0/(nk*sizeof(u4)));
        total_nbytes += nc;
        total_nk += nk;
    }
    printf("\n");
    info("In total, stored %zd kmers in %.2f GB (compression ratio %.2f%%)",
         total_nk, total_nbytes * 1.0/(1024*1024*1024), total_nbytes * 100.0/(total_nk * sizeof(u4)));

    close(fd);
}

