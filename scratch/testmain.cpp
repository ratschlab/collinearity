//
// Created by Sayan Goswami on 29.11.2024.
//

#include "tests.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "../src/utils.h"
#include "../src/index.h"
#include "../src/config.h"
#include "sdsl/vectors.hpp"

struct rf_config_t : config_t {
    int &n_threads = kwarg("num-threads", "number of threads").set_default(1);
    static rf_config_t init(int argc, char *argv[]) {
        auto c = argparse::parse<rf_config_t>(argc, argv);
        c.print();
        return c;
    }
};

int main(int argc, char *argv[]) {
    fna4();
}

fn(a0) {
    string seq =
            "TGTAACCTCCATGTGATGATCTAAAACAATAACAAATAAATAGTTCCTCCCATATAATAT"
            "TATTTCTTACATAATAAAGAATATCATATATTCTCAAAAAATAACAAATAATATCCTCTT"
            "TCCATTCTCAATTAAGTTCTTAAATGAGAATAAAAGGGTAATCCTCTGTATTTCTTAA";

    auto seq1 = create_kmers_1t(seq, 10, 4, encode_dna);
    auto seq2 = create_kmers(seq, 10, 4, encode_dna);
    _verify(parlay::equal(seq1, seq2));
    prettyPrintVector(seq1);
    prettyPrintVector(seq2);
}

fn(a1) {
    auto data = parlay::tabulate(1000, [](size_t i) {return i;});
    auto slice = parlay::make_slice(data.begin() + 10, data.begin() + 30);
    cout << slice.begin() << "\n";
    cout << slice.end() << "\n";
    cout << slice.size() << "\n";
    cout << slice[5] << "\n";
}

fn(a2) {
    string seq =
            "TGTAACCTCCATGTGATGATCTAAAACAATAACAAATAAATAGTTCCTCCCATATAATAT"
            "TATTTCTTACATAATAAAGAATATCATATATTCTCAAAAAATAACAAATAATATCCTCTT"
            "TCCATTCTCAATTAAGTTCTTAAATGAGAATAAAAGGGTAATCCTCTGTATTTCTTAA";
    auto slice = parlay::make_slice(seq.rbegin(), seq.rend());
    auto rc = parlay::delayed_map(slice, [](char c) {
        return "TGAC"[(c >> 1) & 3];
    });
    cout << slice.size() << "\n";
    cout << slice[0] << "\n";
    cout << *slice.begin() << "\n";

}

fn(a3) {
    sdsl::int_vector<> v = {3,2,1,0,2,1,3,4,1,1,1,3,2,3};
    v[1] = 0;
    sdsl::util::bit_compress(v);
    log_info("V:");
    cout << v << endl;
    log_info("Size:");
    cout << sdsl::size_in_bytes(v) << endl;
}

fn(a4) {
    auto filename = "data.bin";
    int a1=10, b1=21, c1=32, a2, b2, c2;
    float d1=0.23f, e1=3.14f, f1=1000.0f, d2, e2, f2;
    bool g1= true, h1= false, i1= true, g2, h2, i2;

    auto fp = fopen(filename, "w");
    dump_values(fp, a1, b1, c1, d1, e1, f1, g1, h1, i1);
    fclose(fp);

    fp = fopen(filename, "r");
    load_values(fp, &a2, &b2, &c2, &d2, &e2, &f2, &g2, &h2, &i2);
    log_info("%d, %d, %d, %.2f, %.2f, %.2f, %d, %d, %d",
             a2, b2, c2, d2, e2, f2, g2, h2, i2);
}