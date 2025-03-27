//
// Created by Sayan Goswami on 29.11.2024.
//

#include "tests.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "../src/utils.h"
#include "../src/index.h"

using namespace std;

int main(int argc, char *argv[]) {
    fa3();
}

f(a0) {
    string seq =
            "TGTAACCTCCATGTGATGATCTAAAACAATAACAAATAAATAGTTCCTCCCATATAATAT"
            "TATTTCTTACATAATAAAGAATATCATATATTCTCAAAAAATAACAAATAATATCCTCTT"
            "TCCATTCTCAATTAAGTTCTTAAATGAGAATAAAAGGGTAATCCTCTGTATTTCTTAA";

    auto seq1 = create_kmers_1t(seq, 10, 4, encode_dna);
    auto seq2 = create_kmers(seq, 10, 4, encode_dna);
    verify(parlay::equal(seq1, seq2));
    prettyPrintVector(seq1);
    prettyPrintVector(seq2);
}

f(a1) {
    auto data = parlay::tabulate(1000, [](size_t i) {return i;});
    auto slice = parlay::make_slice(data.begin() + 10, data.begin() + 30);
    cout << slice.begin() << "\n";
    cout << slice.end() << "\n";
    cout << slice.size() << "\n";
    cout << slice[5] << "\n";
}

f(a2) {
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


void foo(index_t *idx) {
    idx->init_query_buffers();
}