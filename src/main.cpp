#include "collinearity.h"

int main(int argc, char *argv[]) {
    config_t config(argc, argv);

    index_t *idx;
    if (config.jaccard) idx = new j_index_t(config.k, config.sigma, config.fwd_rev);
    else idx = new c_index_t(config.k, config.sigma, config.fwd_rev);

    if (config.phase == config_t::index) {
        if (!config.raw) {
            index_fasta(config.ref.c_str(), idx);
            idx->dump(config.idx);
        } else throw "Not implemented";
    } else if (config.phase == config_t::query) {
        idx->load(config.idx);
        query_fasta(idx, config.qry.c_str(), 4096);
    } else if (config.phase == config_t::both) {
        index_fasta(config.ref.c_str(), idx);
        query_fasta(idx, config.qry.c_str(), 4096);
    }

    return 0;
}