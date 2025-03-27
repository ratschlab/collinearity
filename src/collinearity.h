//
// Created by Sayan Goswami on 27.11.2024.
//

#ifndef COLLINEARITY_COLLINEARITY_H
#define COLLINEARITY_COLLINEARITY_H

#include "prelude.h"
#include <vector>
#include "parlay_utils.h"
#include "kseq++/kseq++.hpp"
#include "cqutils.h"
#include "index.h"
#include "rawsignals.h"
#include "config.h"
#include "utils.h"

void index_fasta(const char* fasta_filename, index_t *idx);

//void index_fasta_raw(const char* fasta_filename, std::string poremodel, index_t *idx);

void query_fasta(index_t *idx, const char* fasta_filename, int batch_sz);

//void query_blow5(index_t *idx, const char* blow5_filename, int batch_sz);

//std::vector<std::tuple<std::string, bool, u4, float>> query_batch(index_t *idx, std::vector<std::string> &sequences);

/**
 * From fasta, generate a crude index
 * @param fasta_filename name of the fasta file
 * @param k k-mer size
 * @param sigma alphabet size
 */
//sindex_t process_fasta(const char* fasta_filename, int k, int sigma);

//index_t process_fasta_jaccard(const char* fasta_filename, int k, int sigma);

//index_t process_fasta_raw(const char* fasta_filename, int k, int sigma, std::string poremodel);

/**
 * query fasta input in index
 * @param filename fasta file
 * @param k kmer length
 * @param sigma alphabet size
 * @param batch_sz batch size
 * @param index index
 * @param refnames reference headers
 */
//void query(const char *filename, int k, int sigma, size_t batch_sz, index_t &index);

//void query_jaccard(const char *filename, int k, int sigma, size_t batch_sz, index_t &index);

//void query_raw(const char *filename, int k, int sigma, size_t batch_sz, index_t &index);


#endif //COLLINEARITY_COLLINEARITY_H
