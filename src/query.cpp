//
// Created by Sayan Goswami on 09.12.2024.
//

#include "collinearity.h"
#include "../external/slow5lib/include/slow5/slow5.h"

#define SANITY_CHECKS 1

using namespace klibpp;

const char STRAND[2] = {'-', '+'};

std::vector<std::tuple<std::string, bool, u4, float>> query_batch(index_t *idx, std::vector<std::string> &sequences) {
    auto nr = sequences.size();
    std::vector<std::tuple<std::string, bool, u4, float>> results;
    parlay::for_each(parlay::iota(nr), [&](size_t i){
        auto result = idx->search(sequences[i]);
        return std::make_tuple(
                std::string(std::get<0>(result)), std::get<1>(result), std::get<2>(result), std::get<3>(result));
    });
    return results;
}

void query_fasta(index_t *idx, std::string &fasta_filename, int batch_sz, std::string &outfile) {
    auto fp = fopen(outfile.c_str(), "w");
    if (!fp) log_error("Could not open %s because %s.", outfile.c_str(), strerror(errno));

    idx->init_query_buffers();
    std::vector<std::string> headers, sequences;
    headers.reserve(batch_sz);
    sequences.reserve(batch_sz);

    // read sequences and convert to key-value pairs
    KSeq record;
    auto fd = open(fasta_filename.c_str(), O_RDONLY);
    if (fd < 0) log_error("Could not open %s because %s.", fasta_filename.c_str(), strerror(errno));
    auto ks = make_kstream(fd, read, mode::in);
    log_info("Begin query..");
    int nr = 0;
    u8 total_nr = 0;
    while (ks >> record) {
        headers.emplace_back(record.name);
        sequences.emplace_back(record.seq);
        nr++;
        if (nr == batch_sz) {
            auto results = parlay::tabulate(nr, [&](size_t i) {
                return idx->search(sequences[i]);
            });
            for (u4 i = 0; i < nr; ++i)
                fprintf(fp, "%s\t%zu\t%s\t%c\t%d\t%f\n", headers[i].c_str(), sequences[i].length(), std::get<0>(results[i]),
                       STRAND[(int)std::get<1>(results[i])], std::get<2>(results[i]), std::get<3>(results[i]));
            total_nr += nr;
            sitrep("%lu", total_nr);
            nr = 0;
            headers.clear();
            sequences.clear();
        }
    }
    if (nr) {
        auto results = parlay::tabulate(nr, [&](size_t i) {
            return idx->search(sequences[i]);
        });
        for (u4 i = 0; i < nr; ++i)
            fprintf(fp, "%s\t%zu\t%s\t%c\t%d\t%f\n", headers[i].c_str(), sequences[i].length(), std::get<0>(results[i]),
                   STRAND[(int)std::get<1>(results[i])], std::get<2>(results[i]), std::get<3>(results[i]));
        total_nr += nr;
        sitrep("%lu", total_nr);
        nr = 0;
        headers.clear();
        sequences.clear();
    }
    stderrflush;
    close(fd);
    fclose(fp);
    log_info("Done.");
}

void query_blow5(index_t *idx, const char* blow5_filename, int batch_sz) {
    throw "Not implemented";
}

