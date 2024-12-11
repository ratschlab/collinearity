#include "collinearity.h"

int main(int argc, char *argv[]) {
    const char *ref = "/scratch/Zymo/Refs-fwd-rev.fasta";
    const char *qry = "/scratch/Zymo/basecalled/fast/Sigs.fasta"; //"/scratch/Zymo/reads-tiny.fasta";
    auto [idx, refnames] = process_fasta(argv[1], 10, 4);
    query(argv[2], 10, 4, 4096, idx, refnames);
    return 0;
}
