//
// Created by Sayan Goswami on 27.11.2024.
//

#include <vector>
#include "prelude.h"
#include <parlay/primitives.h>
#include <parlay/io.h>
#include "kseq++/kseq++.hpp"

using namespace klibpp;

#define encode_dna(x) (((x)>>1)&3)

parlay::sequence<KeyT> create_kmers(const std::string& sequence, int k, int sigma) {
    size_t n = sequence.size();
    if (k > n) return {};
    auto indices = parlay::iota(n - k + 1);
    auto kmers = parlay::map(indices, [&](size_t i) {
        KeyT kmer = 0;
        for (int j = 0; j < k; ++j) kmer = kmer * sigma + encode_dna(sequence[i+j]);
        return kmer;
    });
    return kmers;
}

void process_fasta(const char* fasta_filename, int k, int sigma) {
    // read sequences and convert to key-value pairs
    int ref_id = 0;
    std::vector<std::string> header_names;
    std::vector<uint> ref_lens;
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    auto ks = make_kstream(fd, read, mode::in);

    while (ks >> record) {
        header_names.emplace_back(record.name);
        auto kmers = create_kmers(record.seq, k, sigma);
        auto kmer_str = parlay::to_chars(kmers);
        info("%s", record.name.c_str());
        std::cout << kmer_str << std::endl;
    }

    close(fd);
}

