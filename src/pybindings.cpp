//
// Created by Sayan Goswami on 16.01.2025.
//

#include "collinearity.h"
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

struct alignment_t {
    std::string trg_name = "*";
    int trg_pos = -1;
    float presence = 0.0f;
};

struct rt_index_t {
    int k;
    index_t index;

    explicit rt_index_t(int k): k(k), index(k) {}

    rt_index_t(const std::string& filename, int k):
        k(k), index(process_fasta(filename.c_str(), k, 4)) {}

    rt_index_t(const std::string& filename, int k, std::string &poremodel):
        k(k), index(process_fasta_raw(filename.c_str(), k, 16, poremodel)) {}

    void load(const std::string& basename) {
        index.load(basename);
    }

    void dump(const std::string& basename) {
        index.dump(basename);
    }

    std::vector<alignment_t> query_batch(std::vector<std::string> &sequences) {
        const auto B = sequences.size();
        std::vector<alignment_t> alignments(B);
        parlay::for_each(parlay::iota(B), [&](size_t i){
            const auto &seq = sequences[i];
            u4 n = seq.size();
            parlay::sequence<u4> keys;
            for (int j = 0; j < n-k+1; ++j) {
                keys[j] = encode_kmer(seq.data() + i, k, 4, encode_dna);
                const auto [trg_id, trg_pos, support] = index.search(keys.data(), keys.size());
                if (trg_id == -1) {
                    alignments[i].trg_name = '*', alignments[i].trg_pos = -1, alignments[i].presence = 0;
                } else {
                    alignments[i].trg_name = index.headers[trg_id], alignments[i].trg_pos = trg_pos, alignments[i].presence = support;
                }
            }
        });
        return alignments;
    }

    std::vector<alignment_t> query_batch(std::vector<std::vector<double>> &signals) {
        const auto B = signals.size();
        std::vector<alignment_t> alignments(B);
        parlay::for_each(parlay::iota(B), [&](size_t i){
            auto &signal = signals[i];
            parlay::sequence<double> calibrated_signal(signal.begin(), signal.end());
            tstat_segmenter_t segmenter;
            auto events = generate_events(calibrated_signal, segmenter);
            auto quantized = quantize_signal(events);
            parlay::sequence<u4> keys(quantized.size() - k + 1);
            for (int j = 0; j < quantized.size() - k + 1; ++j)
                keys[i] = encode_kmer(quantized.data() + j, k, 16, encode_qsig);
            const auto [trg_id, trg_pos, support] = index.search(keys.data(), keys.size());
            if (trg_id == -1) {
                alignments[i].trg_name = '*', alignments[i].trg_pos = -1, alignments[i].presence = 0;
            } else {
                alignments[i].trg_name = index.headers[trg_id], alignments[i].trg_pos = trg_pos, alignments[i].presence = support;
            }
        });
        return alignments;
    }

    std::vector<alignment_t> query_batch(std::vector<std::vector<int>> &signals, std::vector<double> &digitisations,
                                         std::vector<double> &offsets, std::vector<double> &ranges) {
        const auto B = signals.size();
        std::vector<alignment_t> alignments(B);
        parlay::for_each(parlay::iota(B), [&](size_t i){
            auto &signal = signals[i];
            parlay::sequence<double> calibrated_signal(signal.size());
            for (int j = 0; j < signal.size(); ++j) {
                calibrated_signal[j] = TO_PICOAMPS(signal[j], digitisations[i], offsets[i], ranges[i]);
            }
            tstat_segmenter_t segmenter;
            auto events = generate_events(calibrated_signal, segmenter);
            auto quantized = quantize_signal(events);
            parlay::sequence<u4> keys(quantized.size() - k + 1);
            for (int j = 0; j < quantized.size() - k + 1; ++j)
                keys[i] = encode_kmer(quantized.data() + j, k, 16, encode_qsig);
            const auto [trg_id, trg_pos, support] = index.search(keys.data(), keys.size());
            if (trg_id == -1) {
                alignments[i].trg_name = '*', alignments[i].trg_pos = -1, alignments[i].presence = 0;
            } else {
                alignments[i].trg_name = index.headers[trg_id], alignments[i].trg_pos = trg_pos, alignments[i].presence = support;
            }
        });
        return alignments;
    }


};

// Binding the function to the Python module
PYBIND11_MODULE(_core, m) {

    pybind11::class_<rt_index_t>(m, "rt_index_t")
        .def(py::init<int>())
        .def(py::init<std::string&, int>())
        .def(py::init<std::string&, int, std::string&>())
        .def("load", &rt_index_t::load)
        .def("dump", &rt_index_t::dump)
        .def("query_batch", static_cast<std::vector<alignment_t> (rt_index_t::*)(std::vector<std::string>&)>(&rt_index_t::query_batch))
        .def("query_batch", static_cast<std::vector<alignment_t> (rt_index_t::*)(std::vector<std::vector<double>>&)>(&rt_index_t::query_batch))
        .def("query_batch", static_cast<std::vector<alignment_t> (rt_index_t::*)(std::vector<std::vector<int>>&, std::vector<double>&,
                                                                                 std::vector<double>&, std::vector<double>&)>(&rt_index_t::query_batch));

    pybind11::class_<alignment_t>(m, "alignment_t")
            .def_readonly("trg_name", &alignment_t::trg_name)
            .def_readonly("trg_pos", &alignment_t::trg_pos)
            .def_readonly("presence", &alignment_t::presence);
}