//
// Created by Sayan Goswami on 27.01.2025.
//

#ifndef COLLINEARITY_CONFIG_H
#define COLLINEARITY_CONFIG_H

#include "prelude.h"
#include "argparse/argparse.hpp"

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

template <class T>
static void dump_values(FILE *fp, T& var) {
    fwrite(&var, sizeof(var), 1, fp);
}

template <class T, typename... Args>
static void dump_values(FILE *fp, T& var, Args... args) {
    expect(fwrite(&var, sizeof(var), 1, fp) == 1);
    dump_values(fp, args...);
}

template <class T>
static void load_values(FILE *fp, T *p_var) {
    expect(fread(p_var, sizeof(*p_var), 1, fp) == 1);
}

template <class T, typename... Args>
static void load_values(FILE *fp, T *p_var, Args... args) {
    expect(fread(p_var, sizeof(*p_var), 1, fp) == 1);
    load_values(fp, args...);
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
    bool &compressed = flag("compressed", "Use a compressed jaccard index.");
    bool &fwd_rev = flag("fr", "Index both forward and reverse strands of the reference.");
    float &presence_fraction = kwarg("pf", "Fraction of k-mers that must be present in an alignment.").set_default(0.1f);
    std::string &_sort_block_size = kwarg("sort-blksz", "Block size to use in sorting.").set_default("");
    int &bandwidth = kwarg("bw", "Width of the band in which kmers contained will be considered collinear").set_default(15);
    int &jc_frag_len = kwarg("jc-frag-len", "If --jaccard is set, the sequence are indexed and queried in overlapping fragments of this length.").set_default(180);
    int &jc_frag_ovlp_len = kwarg("jc-frag-ovlp-len", "If --jaccard is set, the sequence are indexed and queried in fragments which overlap this much.").set_default(120);
    bool &dynamic = flag("dynamic", "use a dynamic multi-map");
    int &n_shard_bits = kwarg("num-shard-bits", "log2(x), where x is the number of shards").set_default(10);
    std::string &poremodel = kwarg("poremodel", "Pore-model file path (currently unused)").set_default("");

    bool raw = false;
    int sigma = 4;
    u8 sort_block_size;

    config_t() = default;

    static config_t init(int argc, char *argv[]) {
        auto c = argparse::parse<config_t>(argc, argv);

        if (c.is_valid) {
            if (c._sort_block_size.empty()) c.sort_block_size = MEMPOOL_BLOCKSZ;
            else c.sort_block_size = hmsize2bytes(c._sort_block_size);
            c.sort_block_size = (c.sort_block_size < MEMPOOL_BLOCKSZ)? MEMPOOL_BLOCKSZ : c.sort_block_size;
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

    void dump(FILE *fp) {
        dump_values(fp, k, jaccard, fwd_rev, presence_fraction, bandwidth, jc_frag_len, jc_frag_ovlp_len, dynamic,
                    n_shard_bits, raw, sigma, sort_block_size);
    }

    void load(FILE *fp) {
        load_values(fp, &k, &jaccard, &fwd_rev, &presence_fraction, &bandwidth, &jc_frag_len, &jc_frag_ovlp_len, &dynamic,
                    &n_shard_bits, &raw, &sigma, &sort_block_size);
    }

private:
    const char *phase_strs[3] = {"index", "query", "both"};
    const char *true_false[2] = {"false", "true"};
};

#endif //COLLINEARITY_CONFIG_H
