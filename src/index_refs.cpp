//
// Created by Sayan Goswami on 27.11.2024.
//

#include "collinearity.h"

#define SANITY_CHECKS 1

using namespace klibpp;

void index_fasta(const char* fasta_filename, index_t *idx) {
    log_info("Beginning indexing..");
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    if (fd < 0) log_error("Could not open %s because %s.", fasta_filename, strerror(errno));
    auto ks = make_kstream(fd, read, mode::in);

    u4 ref_id = 0;
    while (ks >> record) {
        idx->add(record.name, record.seq);
        sitrep("processed %u references.", ++ref_id);
    }

    stderrflush;
    close(fd);
    idx->build();
}

void index_fasta_raw(const char* fasta_filename, std::string poremodel, index_t *idx) {
    throw "Not implemented";
    const auto [pore_k, pore_levels] = load_pore_model(poremodel);
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    if (fd < 0) log_error("Could not open %s because %s.", fasta_filename, strerror(errno));
    auto ks = make_kstream(fd, read, mode::in);

    u4 ref_id = 0;
    while (ks >> record) {
        auto squiggles = sequence2squiggles(record.seq, pore_k, pore_levels);
        auto quant = quantize_signal_simple(squiggles);
        // todo - add to index
        sitrep("processed %u references.", ++ref_id);
    }

    stderrflush;
    close(fd);
    idx->build();
}

