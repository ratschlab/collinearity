#include "collinearity.h"

int main(int argc, char *argv[]) {
    info("Hello!");
    auto config = config_t(argc, argv);
    config.print();

    if (config.phase == 1) {
        // index and dump
        if (config.raw) {
            auto idx = process_fasta_raw(config.ref.c_str(), config.k, config.sigma, config.poremodel);
            idx.calc_max_occ();
            idx.dump(config.ref);
        } else {
            auto idx = process_fasta(config.ref.c_str(), config.k, config.sigma);
            idx.print_info();
            idx.dump(config.ref);
        }
    } else if (config.phase == 2) {
        // load and query
        if (config.raw) {
            sindex_t idx;
            idx.load(config.ref);
            query_raw(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        } else {
            index_t idx(config.k, config.sigma);
            idx.load(config.ref);
            query(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        }
    } else if (config.phase == 3) {
        if (config.raw) {
            auto idx = process_fasta_raw(config.ref.c_str(), config.k, config.sigma, config.poremodel);
            idx.calc_max_occ();
            query_raw(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        } else {
            auto idx = process_fasta(config.ref.c_str(), config.k, config.sigma);
            idx.print_info();
            query(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        }
    }
    info("Bye!");
    return 0;
}