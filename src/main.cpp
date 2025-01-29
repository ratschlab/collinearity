#include "collinearity.h"

int main(int argc, char *argv[]) {
    info("Hello!");
    auto config = config_t(argc, argv);
    config.print();

    if (config.phase == 1) {
        // index and dump
        if (config.raw) {
            auto idx = process_fasta_raw(config.ref.c_str(), config.k, config.sigma, config.poremodel);
            idx.print_info();
            idx.dump(config.ref);
        } else {
            auto idx = process_fasta(config.ref.c_str(), config.k, config.sigma);
            idx.print_info();
            idx.dump(config.ref);
        }
    } else if (config.phase == 2) {
        // load and query
        index_t idx(config.k, config.sigma);
        idx.load(config.ref);
        if (config.raw) {
            query_raw(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        } else {
            query(config.qry.c_str(), config.k, config.sigma, BATCH_SZ, idx);
        }
    }
    info("Bye!");
    return 0;
}