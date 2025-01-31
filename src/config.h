//
// Created by Sayan Goswami on 27.01.2025.
//

#ifndef COLLINEARITY_CONFIG_H
#define COLLINEARITY_CONFIG_H

#include "prelude.h"

using namespace std;

static inline size_t hmsize2bytes(std::string& hsize) {
    char unit = hsize.back();
    if (unit >= 48 && unit <= 57)
        return std::stoul(hsize);
    else {
        hsize.pop_back();
        size_t size = std::stoul(hsize);
        switch (unit) {
            case 'k':
            case 'K':
                size = size KiB;
                break;
            case 'm':
            case 'M':
                size = size MiB;
                break;
            case 'g':
            case 'G':
                size = size GiB;
                break;
            default:
                error("Unknown unit %c.", unit);
        }
        return size;
    }
}

struct config_t {
    int phase = 0;
    string ref, qry, poremodel;
    bool raw = false;
    int k=0, sigma;

    void print() const {
        info("Config:");
        info("Phase: %d", phase);
        info("Ref: %s", ref.c_str());
        info("Query: %s", qry.c_str());
        info("Poremodel: %s, Raw = %d", poremodel.c_str(), raw);
        info("k = %d, sigma = %d", k, sigma);
    }

    const char *usage = R""""(
    collinearity {index/query} [options]

    collinearity index  [fasta]              # path input fasta file \
                        --poremodel []       # if specified, the fasta will be indexed raw using the poremodel file
                        -k                   # kmer length

    collinearity query  [fasta]     # path to fasta file \
                        [query]     # input fasta/fastq/blow5 file. For blow5 file, the ref must be indexed in raw formal
                        -k                   # kmer length (must be the same one used in indexing)

    **Exmaples:**

    collinearity index references.fa -k 15    # creates a `references.fa.cidx` file \
    && collinearity query references.fa reads.fastq -k 15 > report.tsv

    collinearity index references.fa --poremodel models.csv -k 8    # creates a `references.fa.cidx` file \
    && collinearity query references.fa signals.blow5 -k 8 > report.tsv
    )"""";

    void print_usage_and_exit() const {
        printf("%s\n", usage);
        exit(1);
    }

    config_t(const config_t&) = delete;
    config_t& operator = (const config_t&) = delete;

    config_t(int argc, char *argv[]) {
        if (argc >= 3) {
            --argc, ++argv;
            if (streq(argv[0], "index")) {
                phase = 1;
                ref = argv[1];
                if (ref.substr(0,2) == "--") print_usage_and_exit();
                argc -= 2, argv += 2;
            } else if (streq(argv[0], "query")) {
                phase = 2;
                if (argc < 3) print_usage_and_exit();
                ref = argv[1];
                if (ref.substr(0,2) == "--") print_usage_and_exit();
                qry = argv[2];
                if (qry.substr(0,2) == "--") print_usage_and_exit();
                argc -= 3, argv += 3;
                if (str_endswith(qry.c_str(), ".blow5")) raw = true;
            } else if (streq(argv[0], "both")) {
                phase = 3;
                if (argc < 3) print_usage_and_exit();
                ref = argv[1];
                if (ref.substr(0,2) == "--") print_usage_and_exit();
                qry = argv[2];
                if (qry.substr(0,2) == "--") print_usage_and_exit();
                argc -= 3, argv += 3;
            }
            else print_usage_and_exit();
        } else print_usage_and_exit();

        while (argc) {
            if (streq(*argv, "--poremodel")) {
                raw = true;
                poremodel = argv[1];
                if (poremodel.substr(0,2) == "--") print_usage_and_exit();
                argc -= 2, argv += 2;
            } else if (streq(argv[0], "-k")) {
                k = atoi(reinterpret_cast<const char *>(argv[1]));
                argc -= 2, argv += 2;
            } else print_usage_and_exit();
        }

        if (str_endswith(qry.c_str(), ".blow5")) raw = true;

        if (raw) {
            if (k < 1 || k > 10) k = 8;
            sigma = 16;
        }
        else {
            if (k < 1 || k > 16) k = 15;
            sigma = 4;
        }
    }

};

#endif //COLLINEARITY_CONFIG_H
