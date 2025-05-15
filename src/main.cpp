#include "collinearity.h"

int main(int argc, char *argv[]) {
    auto config = config_t::init(argc, argv);

    index_t *idx;
    if (config.phase == config_t::index) {
        if (config.jaccard) idx = new j_index_t(config);
        else idx = new c_index_t(config);
        if (!config.raw) {
            index_fasta(config.ref, idx);
            dump_index(config.idx, config, idx);
        } else throw "Not implemented";
    } else if (config.phase == config_t::query) {
        idx = load_index(config.idx);
        query_fasta(idx, config.qry, 4096, config.out);
    } else if (config.phase == config_t::both) {
        if (config.jaccard) {
            if (config.compressed) idx = new cj_index_t(config);
            else idx = new j_index_t(config);
        }
        else idx = new c_index_t(config);
        index_fasta(config.ref, idx);
        query_fasta(idx, config.qry, 4096, config.out);
    }

    return 0;
}