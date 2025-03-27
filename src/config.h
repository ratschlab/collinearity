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

struct config_t : public argparse::Args {
    enum phase_t {index, query, both};
    phase_t phase;
    std::string idx = "";

    std::string &ref = kwarg("ref", "Path of input reference fasta file. If this is not set, then the --idx must be set.").set_default("");
    std::string &_idx = kwarg("idx", "Base path to output index file (a .cidx suffix will be added to the path). "
                                    "If --idx and --qry are not provided, the --ref path is used to dump the index."
                                    "If this is not set, then the --ref file must be set."
                                    "This is ignored if --ref path is provided.").set_default("");
    std::string &qry = kwarg("qry", "Path to query fasta. If not provided, then the index is build and dumped to file.").set_default("");
    std::string &out = kwarg("out", "Path to output file. Must be provided if --qry is provided.").set_default("");
    int &k = kwarg("k", "k-mer length").set_default(15);
    bool &jaccard = flag("jaccard", "Use jaccard similarity.");
    bool &fwd_rev = flag("fr", "Index both forward and reverse strands of the reference.");
    float &presence_fraction = kwarg("pf", "Fraction of k-mers that must be present in an alignment.").set_default(0.3f);
    int &bandwidth = kwarg("bw", "Width of the band in which kmers contained will be considered collinear").set_default(15);
    int &jc_frag_len = kwarg("jc-frag-len", "If --jaccard is set, the sequence are indexed and queried in overlapping fragments of this length.").set_default(180);
    int &jc_frag_ovlp_len = kwarg("jc-frag-ovlp-len", "If --jaccard is set, the sequence are indexed and queried in fragments which overlap this much.").set_default(120);
    std::string &poremodel = kwarg("poremodel", "Pore-model file path (currently unused)").set_default("");

    bool raw = false;
    int sigma = 4;

    config_t() = default;

    static config_t init(int argc, char *argv[]) {
        auto c = argparse::parse<config_t>(argc, argv);

        if (c.is_valid) {
            c.idx = c._idx + ".cidx";
            if (c.ref.empty() && c._idx.empty()) goto PRINT_HELP_AND_EXIT;
            if (c.ref.empty() && c.qry.empty()) goto PRINT_HELP_AND_EXIT;
            if (c._idx.empty() && c.qry.empty()) {
                c.idx = c.ref + ".cidx";
                c.phase = index;
            } else if (c.qry.empty()) {
                c.phase = index;
            } else if (!c.ref.empty() && !c.qry.empty()) {
                c.phase = both;
                c.idx = "";
                if (c.out.empty()) goto PRINT_HELP_AND_EXIT;
            } else if (!c._idx.empty() && !c.qry.empty()) {
                c.phase = query;
                if (c.out.empty()) goto PRINT_HELP_AND_EXIT;
            }
        } else {
            PRINT_HELP_AND_EXIT:
            c.help();
            exit(1);
        }

        c.print();
        fprintf(stderr, "Phase: %s\n", c.phase_strs[c.phase]);
        fprintf(stderr, "Index: %s\n", c.idx.c_str());
        return c;
    }

private:
    const char *phase_strs[3] = {"index", "query", "both"};
    const char *true_false[2] = {"false", "true"};

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
