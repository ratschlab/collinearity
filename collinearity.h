//
// Created by Sayan Goswami on 27.11.2024.
//

#ifndef COLLINEARITY_COLLINEARITY_H
#define COLLINEARITY_COLLINEARITY_H

#include "prelude.h"

/**
 * From fasta, generate a crude index
 * @param fasta_filename name of the fasta file
 * @param k k-mer size
 * @param sigma alphabet size
 */
void * process_fasta(const char* fasta_filename, int k, int sigma);

#endif //COLLINEARITY_COLLINEARITY_H
