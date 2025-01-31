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
#include "rawsignals.h"
#include "config.h"
#include "utils.h"

/**
 * From fasta, generate a crude index
 * @param fasta_filename name of the fasta file
 * @param k k-mer size
 * @param sigma alphabet size
 */
index_t process_fasta(const char* fasta_filename, int k, int sigma);

sindex_t process_fasta_raw(const char* fasta_filename, int k, int sigma, std::string poremodel);

/**
 * query fasta input in index
 * @param filename fasta file
 * @param k kmer length
 * @param sigma alphabet size
 * @param batch_sz batch size
 * @param index index
 * @param refnames reference headers
 */
void query(const char *filename, int k, int sigma, size_t batch_sz, index_t &index);

void query_raw(const char *filename, int k, int sigma, size_t batch_sz, sindex_t &index);

static auto encode_dna = [](char x) { return (u4)(((x) >> 1) & 3); };
static auto encode_qsig = [](u1 x) { return x; };


#endif //COLLINEARITY_COLLINEARITY_H
