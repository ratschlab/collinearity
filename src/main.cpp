#include "collinearity.h"

int main(int argc, char *argv[]) {
    info("Hello!");
    const char *ref = "/scratch/Zymo/Refs-fwd-rev.fasta";
    const char *qry = "/scratch/Zymo/basecalled/fast/Sigs.fasta"; //"/scratch/Zymo/reads-tiny.fasta";

    if (str_endswith(argv[1], ".cidx")) {
        auto [idx, refnames] = load_index(argv[1]);
        if (SANITY_CHECKS) print_idx_info(idx);
        query(argv[2], KMER_LENGTH, SIGMA, BATCH_SZ, idx, refnames);
    } else {
        auto [idx, refnames] = process_fasta(argv[1], KMER_LENGTH, SIGMA);
        dump_index(idx, refnames, argv[1]);
        if (SANITY_CHECKS) print_idx_info(idx);
        if (argc > 2) query(argv[2], KMER_LENGTH, SIGMA, BATCH_SZ, idx, refnames);
    }
    info("Bye!");
    return 0;
}