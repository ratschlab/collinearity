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

void index_fasta(std::string &fasta_filename, index_t *idx);

void query_fasta(index_t *idx, std::string &asta_filename, int batch_sz, std::string &outfile);


#endif //COLLINEARITY_COLLINEARITY_H
