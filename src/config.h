//
// Created by Sayan Goswami on 27.01.2025.
//

#ifndef COLLINEARITY_CONFIG_H
#define COLLINEARITY_CONFIG_H

#include "prelude.h"
#include "argparse/argparse.hpp"

using namespace std;

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

struct config_t : public argparse::Args {
    enum phase_t {index, query, both};
    struct idx_args_t : public argparse::Args {
        std::string &ref = kwarg("ref", "reference file path");
        std::string &poremodel = kwarg("poremodel", "pore-model file path").set_default("");
        int &k = kwarg("k", "k-mer length").set_default(15);
    };

    struct qry_args_t : public argparse::Args {
        std::string &ref = kwarg("ref", "reference file path");
        std::string &qry = kwarg("qry", "query file path");
    };

    struct idx_qry_args_t : public argparse::Args {
        std::string &ref = kwarg("ref", "reference file path");
        std::string &poremodel = kwarg("poremodel", "pore-model file path").set_default("");
        int &k = kwarg("k", "k-mer length").set_default(15);
        std::string &qry = kwarg("qry", "query file path");
    };

    idx_args_t &idx_args = subcommand("index");
    qry_args_t &qry_args = subcommand("query");
    idx_qry_args_t &both_args = subcommand("both");

    std::string ref, qry, poremodel;
    bool raw = false;
    int k = 0, sigma = 0;
    phase_t phase = both;

    config_t init(int argc, char *argv[]) {
        auto config = argparse::parse<config_t>(argc, argv);
        if (config.idx_args.is_valid) {
            ref = config.idx_args.ref;
            poremodel = config.idx_args.poremodel;
            k = config.idx_args.k;
            raw = !poremodel.empty();
        } else if (config.qry_args.is_valid) {
            ref = config.qry_args.ref;
            qry = config.qry_args.qry;
            raw = str_endswith(qry.c_str(), ".blow5");
        } else {
            ref = config.both_args.ref;
            poremodel = config.both_args.poremodel;
            k = config.both_args.k;
            qry = config.both_args.qry;
            raw = !poremodel.empty();
        }
        if (raw) {
            if (k < 1 || k > 10) k = 8;
            sigma = 16;
        }
        else {
            if (k < 1 || k > 16) k = 15;
            sigma = 4;
        }
        return config;
    }
};

#endif //COLLINEARITY_CONFIG_H
