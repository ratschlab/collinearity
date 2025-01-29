//
// Created by Sayan Goswami on 27.11.2024.
//

#include "collinearity.h"

#define SANITY_CHECKS 1

using namespace klibpp;

index_t process_fasta(const char* fasta_filename, int k, int sigma) {
    index_t index(k, sigma);
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    if (fd < 0) error("Could not open %s because %s.", fasta_filename, strerror(errno));
    auto ks = make_kstream(fd, read, mode::in);

    size_t total_nk = 0, total_nbytes = 0;
    uint64_t ref_id = 0;
    while (ks >> record) {
        if (record.seq.size() > k) {
            auto kmers = create_kmers(record.seq, k, sigma, encode_dna);
            const size_t nk = kmers.size();
            auto addresses = create_addresses(ref_id, nk);
            std::string name = record.name + "+";
            index.add(name, kmers, addresses);
            total_nk += nk;
            ref_id++;

            kmers = create_kmers(revcmp_par(record.seq), k, sigma, encode_dna);
            expect(nk == kmers.size());
            addresses = create_addresses(ref_id, nk);
            name = record.name + "-";
            index.add(name, kmers, addresses);
            total_nk += nk;
            ref_id++;

            sitrep("Generated %zd tuples from %zd sequences", total_nk, ref_id);
        }
    }
    stderrflush;
    close(fd);
    info("Generated %zd tuples from %zd sequences", total_nk, ref_id);

    index.build();
    return index;
}

index_t process_fasta_raw(const char* fasta_filename, int k, int sigma, std::string poremodel) {
    const auto [pore_k, pore_levels] = load_pore_model(poremodel);
    index_t index(k, sigma);
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    if (fd < 0) error("Could not open %s because %s.", fasta_filename, strerror(errno));
    auto ks = make_kstream(fd, read, mode::in);

    size_t total_nk = 0, total_nbytes = 0;
    uint64_t ref_id = 0;
    while (ks >> record) {
        if (record.seq.size() > k + pore_k) {
            auto squiggles = sequence2squiggles(record.seq, pore_k, pore_levels);
            auto quant = quantize_signal(squiggles);
            auto kmers = create_kmers(quant, k, sigma, encode_qsig);
            auto nk = kmers.size();
            auto addresses = create_addresses(ref_id, nk);
            std::string name = record.name + "+";
            index.add(name, kmers, addresses);
            total_nk += nk;
            ref_id++;

            squiggles = sequence2squiggles(revcmp_par(record.seq), pore_k, pore_levels);
            quant = quantize_signal(squiggles);
            kmers = create_kmers(quant, k, sigma, encode_qsig);
            nk = kmers.size();
            addresses = create_addresses(ref_id, nk);
            name = record.name + "-";
            index.add(name, kmers, addresses);
            total_nk += nk;
            ref_id++;

            sitrep("Generated %zd tuples from %zd sequences", total_nk, ref_id);
        }
    }
    stderrflush;
    close(fd);
    info("Generated %zd tuples from %zd sequences", total_nk, ref_id);

    index.build();
    return index;
}

