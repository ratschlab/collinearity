#include "collinearity.h"

int main(int argc, char *argv[]) {
    const char *ref = "/scratch/Zymo/Refs-fwd-rev.fasta";
    const char *qry = "/scratch/Zymo/basecalled/fast/Sigs.fasta"; //"/scratch/Zymo/reads-tiny.fasta";
    auto [idx, refnames] = process_fasta(ref, 10, 4);
    query(qry, 10, 4, 512, idx, refnames);
    return 0;
}
