//
// Created by Sayan Goswami on 27.11.2024.
//

#ifndef COLLINEARITY_COLLINEARITY_H
#define COLLINEARITY_COLLINEARITY_H

#include "prelude.h"
#include <vector>
#include "compressed_array.h"
#include "parlay_utils.h"
#include "kseq++/kseq++.hpp"
#include "cqutils.h"
#include "parlay_hash/unordered_map.h"
#include "index.h"

#define encode_dna(x) (((x)>>1)&3)
// I am assuming that the number of sequences won't exceed 2^20 (1M)
// and the longest sequence won't be longer than 2^40 (1T)
// if that's not the case, change the following line accordingly
#define ref_id_nbits 20
#define ref_len_nbits (64 - ref_id_nbits)
#define ref_id_bitmask (((1ULL) << ref_len_nbits) - 1)
#define make_key_from(id, pos) (((u8)(id)) << (ref_len_nbits) | (pos))
#define get_id_from(key) ((key) >> ref_len_nbits)
#define get_pos_from(key) ((key) & ref_id_bitmask)

// the minimum fraction of kmers required for majority vote
const float presence_fraction = .1f;

// this is the bandwidth, I will later set this from the config
const int bandwidth = 15;

/**
 * From fasta, generate a crude index
 * @param fasta_filename name of the fasta file
 * @param k k-mer size
 * @param sigma alphabet size
 */
std::pair<index_t*, std::vector<std::string>> process_fasta(const char* fasta_filename, int k, int sigma);

/**
 * query fasta input in index
 * @param filename fasta file
 * @param k kmer length
 * @param sigma alphabet size
 * @param batch_sz batch size
 * @param index index
 * @param refnames reference headers
 */
void query(const char *filename, int k, int sigma, size_t batch_sz, index_t *index, std::vector<std::string> &refnames);

static parlay::sequence<u4> create_kmers(const std::string& sequence, int k, int sigma) {
    size_t n = sequence.size();
    if (k > n) return {};
    auto indices = parlay::iota(n - k + 1);
    return parlay::map(indices, [&](size_t i) {
        u4 kmer = 0;
        for (int j = 0; j < k; ++j) kmer = kmer * sigma + encode_dna(sequence[i+j]);
        return kmer;
    });
}

static void print_idx_info(index_t *idx) {
    auto counts_copy = parlay::integer_sort(idx->counts);
    auto start = lower_bound<u4>(counts_copy.data(), 0, counts_copy.size(), 1);
    size_t nk = counts_copy.size() - start;
    info("# kmers = %zd", nk);
    info("Min occ. = %u", counts_copy[start]? counts_copy[start] : counts_copy[start+1]);
    info("Median occ. = %u", counts_copy[start + nk / 2]);
    info("Max occ. = %u", counts_copy.back());
}

#endif //COLLINEARITY_COLLINEARITY_H
