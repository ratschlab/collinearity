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

struct rf_config_t : args_t {
    int &n_threads = kwarg("num-threads", "number of threads").set_default(1);
    static rf_config_t init(int argc, char *argv[]) {
        auto c = argparse::parse<rf_config_t>(argc, argv);
        c.print();
        return c;
    }
};

int main(int argc, char *argv[]) {
    fna6();
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

    auto fout = std::ofstream(filename, std::ios::binary);
    dump_values(fout, a1, b1, c1, d1, e1, f1, g1, h1, i1);
    fout.close();

    auto fin = std::ifstream(filename, std::ios::binary);
    load_values(fin, &a2, &b2, &c2, &d2, &e2, &f2, &g2, &h2, &i2);
    log_info("%d, %d, %d, %.2f, %.2f, %.2f, %d, %d, %d",
             a2, b2, c2, d2, e2, f2, g2, h2, i2);
    fin.close();
}

fn(a5) {
    const u4 n = 1<<20;

    // generate parlay sequence with random numbers
    log_info("Generating data..");
    parlay::random_generator gen;
    std::uniform_int_distribution<u4> dis(0, n-1);
    auto data = parlay::tabulate(n, [&](size_t i) {
        auto r = gen[i];
        return dis(r);
    });
    log_info("Original Size (MB) = %.f", n * 4.0 / (1024.0 * 1024.0));

    // copy the numbers to a cqueue
    cqueue_t<u4> cq;
    cq.push_back(data.data(), n);

    // create an enc vector
    log_info("Creating encoded vectors");
    sdsl::enc_vector<> ev1(data);
    sdsl::enc_vector<> ev2(cq);
    log_info("Before sorting..");
    log_info("Size (MB) = %f, %f", sdsl::size_in_mega_bytes(ev1), sdsl::size_in_mega_bytes(ev2));

    // sort data and create enc vectors
    log_info("Sorting data..");
    auto s_data = parlay::sort(data);
    cqueue_t<u4> s_cq;
    s_cq.push_back(s_data.data(), n);

    log_info("Creating encoded vectors");
    ev1 = sdsl::enc_vector<>(s_data);
    ev2 = sdsl::enc_vector<>(s_cq);
    log_info("After sorting..");
    log_info("Size (MB) = %f, %f", sdsl::size_in_mega_bytes(ev1), sdsl::size_in_mega_bytes(ev2));

    // create a vlc vector
    log_info("Creating vlc vectors..");
    sdsl::vlc_vector<> vv1(data);
    sdsl::vlc_vector<> vv2(cq);
    log_info("Before sorting..");
    log_info("Size (MB) = %f, %f", sdsl::size_in_mega_bytes(vv1), sdsl::size_in_mega_bytes(vv2));

    log_info("Creating vlc vectors..");
    vv1 = sdsl::vlc_vector<>(s_data);
    vv2 = sdsl::vlc_vector<>(s_cq);
    log_info("After sorting..");
    log_info("Size (MB) = %f, %f", sdsl::size_in_mega_bytes(vv1), sdsl::size_in_mega_bytes(vv2));
}

fn(a6) {
    const u4 n = 1<<20;

    // generate parlay sequence with random numbers
    log_info("Generating data..");
    parlay::random_generator gen;
    std::uniform_int_distribution<u4> dis(0, n-1);
    auto data1 = parlay::tabulate(n, [&](size_t i) {
        auto r = gen[i];
        return dis(r);
    });

    auto data2 = parlay::tabulate(n, [&](size_t i) {
        auto r = gen[i];
        return dis(r);
    });

    auto filename = "data.bin";

    // create an enc vector
    log_info("Creating encoded vectors");
    sdsl::enc_vector<> ev1(data1), ev2(data2);
    log_info("Size = %f, %f", sdsl::size_in_mega_bytes(ev1), sdsl::size_in_mega_bytes(ev2));

    auto fout = std::ofstream(filename, std::ios::binary);
    ev1.serialize(fout);
    ev2.serialize(fout);
    dump_values(fout, n);
    fout.close();

    auto fin = std::ifstream(filename, std::ios::binary);
    sdsl::enc_vector<> ev3, ev4;
    sdsl::load(ev3, fin);
    sdsl::load(ev4, fin);
    log_info("Size = %f, %f", sdsl::size_in_mega_bytes(ev3), sdsl::size_in_mega_bytes(ev4));

    _verify(ev1.size() == ev3.size());
    _verify(ev2.size() == ev4.size());

    auto same = parlay::tabulate(n, [&](size_t i) {
        return ev1[i] == ev3[i];
    });
    _verify(parlay::all_of(same, [](bool x) {return x;}));

    same = parlay::tabulate(n, [&](size_t i) {
        return ev2[i] == ev4[i];
    });
    _verify(parlay::all_of(same, [](bool x) {return x;}));

    u4 n0 = 0;
    load_values(fin, &n0);
    _verify(n0 == n);
    fin.close();
}
