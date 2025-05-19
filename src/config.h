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

struct args_t: public argparse::Args {
    std::string &ref = kwarg("ref", "Path of input reference fasta file. If this is not set, then the --idx must be set.").set_default("");
    std::string &idx = kwarg("idx", "Base path to output index file (a .cidx suffix will be added to the path). "
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
    std::string &sort_block_size = kwarg("sort-blksz", "Block size to use in sorting.").set_default("");
    int &bandwidth = kwarg("bw", "Width of the band in which kmers contained will be considered collinear").set_default(15);
    int &jc_frag_len = kwarg("jc-frag-len", "If --jaccard is set, the sequence are indexed and queried in overlapping fragments of this length.").set_default(180);
    int &jc_frag_ovlp_len = kwarg("jc-frag-ovlp-len", "If --jaccard is set, the sequence are indexed and queried in fragments which overlap this much.").set_default(120);
    bool &dynamic = flag("dynamic", "use a dynamic multi-map");
    int &n_shard_bits = kwarg("num-shard-bits", "log2(x), where x is the number of shards").set_default(10);
    int &n_threads = kwarg("n_threads", "Number of threads to use (set <=0 to use all cores)").set_default(0);
};

struct config_t {
    enum phase_t {index, query, both};
    phase_t phase;
    std::string ref, idx, qry, out;
    int sigma=4, k, bandwidth, jc_frag_len, jc_frag_ovlp_len, n_shard_bits, n_threads;
    float presence_fraction;
    bool jaccard, compressed, fwd_rev, dynamic;
    u8 sort_block_size;

    config_t() = default;

    config_t(int argc, char *argv[], bool validate=true) {
        auto args = argparse::parse<args_t>(argc, argv);
        init_from_args(args, validate);
    }

    explicit config_t(std::vector<char*> args, bool validate= false): config_t((int)args.size(), args.data(), validate) {}

    void dump_to(std::ostream &f) {
        dump_values(f, k, bandwidth, jc_frag_len, jc_frag_ovlp_len, n_shard_bits, n_threads, presence_fraction,
                    jaccard, compressed, fwd_rev, dynamic, sort_block_size);
    }

    void load_from(std::istream &f) {
        load_values(f, &k, &bandwidth, &jc_frag_len, &jc_frag_ovlp_len, &n_shard_bits, &n_threads, &presence_fraction,
                    &jaccard, &compressed, &fwd_rev, &dynamic, &sort_block_size);
    }

private:
    void init_from_args(args_t &args, bool validate) {
        ref = args.ref, idx = args.idx, qry = args.qry, out = args.out;
        k=args.k, bandwidth=args.bandwidth, jc_frag_len=args.jc_frag_len, jc_frag_ovlp_len=args.jc_frag_ovlp_len,
        n_shard_bits=args.n_shard_bits, n_threads=args.n_threads;
        presence_fraction=args.presence_fraction;
        jaccard=args.jaccard, compressed=args.compressed, fwd_rev=args.fwd_rev, dynamic=args.dynamic;

        if (args.n_threads > 0) setenv("PARLAY_NUM_THREADS", std::to_string(args.n_threads).c_str(), 1);
        if (args.sort_block_size.empty()) sort_block_size = MEMPOOL_BLOCKSZ;
        else sort_block_size = hmsize2bytes(args.sort_block_size);
        sort_block_size = (sort_block_size < MEMPOOL_BLOCKSZ)? MEMPOOL_BLOCKSZ : sort_block_size;

        if (validate and !is_valid()) {
            args.help(); exit(1);
        }
        args.print();
    }

    bool is_valid() {
        if (!ref.empty() && !qry.empty()) {
            phase = config_t::phase_t::both;
            idx = "";
            if (out.empty()) {
                log_warn("Missing argument: `out`");
                return false;
            }
            return true;
        }
        if (!ref.empty() && qry.empty()) {
            phase = config_t::phase_t::index;
            if (idx.empty()) idx = ref;
            idx = idx + ".cidx";
            return true;
        }
        if (!idx.empty() && !qry.empty()) {
            phase = config_t::phase_t::query;
            idx = idx + ".cidx";
            if (out.empty()) {
                log_warn("Missing argument: `out`");
                return false;
            }
            return true;
        }
        return false;
    }
};

#endif //COLLINEARITY_CONFIG_H
