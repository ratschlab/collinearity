#include "collinearity.h"

int main(int argc, char *argv[]) {
    log_info("Hello!");
    config_t config;
    config.init(argc, argv);

    if (config.phase == config_t::index) {
        // index and dump
        if (config.raw) {
            auto idx = process_fasta_raw(config.ref.c_str(), config.k, config.sigma, config.poremodel);
            idx.dump(config.ref);
        } else {
            auto idx = process_fasta(config.ref.c_str(), config.k, config.sigma);
            idx.dump(config.ref);
        }
    } else if (config.phase == config_t::query) {
        // load and query
        if (config.raw) {
            sindex_t idx;
            idx.load(config.ref);
            query_raw(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        } else {
            sindex_t idx;
            idx.load(config.ref);
            query(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        }
    } else if (config.phase == config_t::both) {
        if (config.raw) {
            auto idx = process_fasta_raw(config.ref.c_str(), config.k, config.sigma, config.poremodel);
            query_raw(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        } else {
            auto idx = process_fasta_jaccard(config.ref.c_str(), config.k, config.sigma);
            query_jaccard(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        }
    } else {
        log_error("Incorrect phase. Options: {load, query, both}");
    }
    log_info("Bye!");
    return 0;
}