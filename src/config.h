//
// Created by Sayan Goswami on 27.01.2025.
//

#ifndef COLLINEARITY_CONFIG_H
#define COLLINEARITY_CONFIG_H

#include "prelude.h"
#include "argparse/argparse.hpp"

#define CONFIG_DUMP_SIZE 3

[[maybe_unused]] static inline size_t hmsize2bytes(std::string& hsize) {
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
                log_error("Unknown unit %c.", unit);
        }
        return size;
    }
}

struct args_t : public argparse::Args {
//    struct idx_args_t : public argparse::Args {
//        std::string &ref = arg("ref", "path of input reference fasta file");
//        std::string &idx = kwarg("idx", "base path to output index file (a .cidx suffix will be added to the path). "
//                                        "If empty, the ref path is used").set_default("");
//        int &k = kwarg("k", "k-mer length").set_default(15);
//        bool &jaccard = flag("jaccard", "use jaccard similarity");
//        bool &fwd_rev = flag("fr", "index both forward and reverse strands of the reference");
//        std::string &poremodel = kwarg("poremodel", "pore-model file path").set_default("");
//    };
//
//    struct qry_args_t : public argparse::Args {
//        std::string &idx = arg("idx", "base path to input index file (a .cidx suffix will be added to the path)");
//        std::string &qry = arg("qry", "query file path");
//    };
//
//    idx_args_t &idx_args = subcommand("index");
//    qry_args_t &qry_args = subcommand("query");

    std::string &ref = kwarg("ref", "path of input reference fasta file. If this is not set, then the --idx must be set.").set_default("");
    std::string &idx = kwarg("idx", "base path to output index file (a .cidx suffix will be added to the path). "
                                    "If --idx and --qry are not provided, the --ref path is used to dump the index."
                                    "If this is not set, then the --ref file must be set."
                                    "This is ignored if --ref path is provided.").set_default("");
    std::string &qry = kwarg("qry", "Path to query fasta. If not provided, then the index is build and dumped to file.").set_default("");
    int &k = kwarg("k", "k-mer length").set_default(15);
    bool &jaccard = flag("jaccard", "Use jaccard similarity.");
    bool &fwd_rev = flag("fr", "Index both forward and reverse strands of the reference.");
    std::string &poremodel = kwarg("poremodel", "Pore-model file path (currently unused)").set_default("");

    static args_t init(int argc, char *argv[]) {
        return argparse::parse<args_t>(argc, argv);
    }
};

struct config_t {
    enum phase_t {index, query, both};
    phase_t phase;
    std::string ref, idx, poremodel, qry;
    bool raw = false, fwd_rev = false, jaccard = false;
    int k, sigma = 4;

    config_t(int argc, char *argv[]) {
        auto args = args_t::init(argc, argv);
        if (args.is_valid) {
            ref = args.ref, idx = args.idx, qry = args.qry;
            k = args.k, jaccard = args.jaccard, fwd_rev = args.fwd_rev;
            if (ref.empty() && idx.empty()) goto PRINT_HELP_AND_EXIT;
            if (ref.empty() && qry.empty()) goto PRINT_HELP_AND_EXIT;
            if (idx.empty() && qry.empty()) {
                idx = ref + ".cidx";
                phase = index;
            } else if (qry.empty()) {
                idx = idx + ".cidx";
                phase = index;
            } else if (!ref.empty() && !qry.empty()) {
                phase = both;
                idx = "";
            } else if (!idx.empty() && !qry.empty()) {
                phase = query;
            }
        } else {
            PRINT_HELP_AND_EXIT:
            args.help();
            exit(1);
        }

        print();
    }

private:
    const char *phase_strs[3] = {"index", "query", "both"};
    const char *true_false[2] = {"false", "true"};
    void print() const {
        fprintf(stderr, "#################\n");
        fprintf(stderr, "Phase: %s\n", phase_strs[phase]);
        fprintf(stderr, "Reference: %s\n", ref.c_str());
        fprintf(stderr, "Index: %s\n", idx.c_str());
        fprintf(stderr, "Poremodel: %s\n", poremodel.c_str());
        fprintf(stderr, "Query: %s\n", qry.c_str());
        fprintf(stderr, "Raw: %s\n", true_false[raw]);
        fprintf(stderr, "Forward & reverse: %s\n", true_false[fwd_rev]);
        fprintf(stderr, "Use Jaccard similarity: %s\n", true_false[jaccard]);
        fprintf(stderr, "k-mer length: %d\n", k);
        fprintf(stderr, "alphabet size: %d\n", sigma);
    }

    void dump(std::string &filename) const {
        u1 buffer[CONFIG_DUMP_SIZE] = {0};
        buffer[0] = ((int)raw) | ((int)fwd_rev << 1) | ((int)jaccard << 2);
        buffer[1] = k, buffer[2] = sigma;
        FILE *fp = fopen(filename.c_str(), "w");
        if (!fp) log_error("Could not open file %s because %s.", filename.c_str(), strerror(errno));
        expect(fwrite(buffer, sizeof(u1), CONFIG_DUMP_SIZE, fp) == CONFIG_DUMP_SIZE);
        fclose(fp);
    }

    void load(std::string &filename) {
        u1 buffer[CONFIG_DUMP_SIZE] = {0};
        FILE *fp = fopen(filename.c_str(), "r");
        if (!fp) log_error("Could not open file %s because %s.", filename.c_str(), strerror(errno));
        expect(fread(buffer, sizeof(u1), CONFIG_DUMP_SIZE, fp) == CONFIG_DUMP_SIZE);
        fclose(fp);
        raw = buffer[0] & 0b1, fwd_rev = buffer[0] & 0b10, jaccard = buffer[0] & 0b100;
        k = buffer[1], sigma = buffer[2];
    }
};

#endif //COLLINEARITY_CONFIG_H
