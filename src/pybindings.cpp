//
// Created by Sayan Goswami on 16.01.2025.
//

#include "collinearity.h"
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

struct Alignment {
    std::string ctg;
    int r_st, r_en, strand;
    float pres_frac;

    Alignment() = default;

    Alignment(const char *header, bool fwd, int start, float pres_frac, int qry_len):
        ctg(header), r_st(start), r_en(start + qry_len), strand(fwd?1:-1), pres_frac(pres_frac) {}
};

//struct rt_index_t {
//    int k;
//    index_t index;
//
//    explicit rt_index_t(int k): k(k), index(k) {}
//
//    rt_index_t(const std::string& filename, int k):
//        k(k), index(process_fasta(filename.c_str(), k, 4)) {}
//
//    rt_index_t(const std::string& filename, int k, std::string &poremodel):
//        k(k), index(process_fasta_raw(filename.c_str(), k, 16, poremodel)) {}
//
//    void load(const std::string& basename) {
//        index.load(basename);
//    }
//
//    void dump(const std::string& basename) {
//        index.dump(basename);
//    }
//
//    std::vector<alignment_t> query_batch(std::vector<std::string> &sequences) {
//        const auto B = sequences.size();
//        std::vector<alignment_t> alignments(B);
//        parlay::for_each(parlay::iota(B), [&](size_t i){
//            const auto &seq = sequences[i];
//            u4 n = seq.size();
//            parlay::sequence<u4> keys;
//            for (int j = 0; j < n-k+1; ++j) {
//                keys[j] = encode_kmer(seq.data() + i, k, 4, encode_dna);
//                const auto [trg_id, trg_pos, support] = index.search(keys.data(), keys.size());
//                if (trg_id == -1) {
//                    alignments[i].trg_name = '*', alignments[i].trg_pos = -1, alignments[i].presence = 0;
//                } else {
//                    alignments[i].trg_name = index.headers[trg_id], alignments[i].trg_pos = trg_pos, alignments[i].presence = support;
//                }
//            }
//        });
//        return alignments;
//    }
//
//    std::vector<alignment_t> query_batch(std::vector<std::vector<double>> &signals) {
//        const auto B = signals.size();
//        std::vector<alignment_t> alignments(B);
//        parlay::for_each(parlay::iota(B), [&](size_t i){
//            auto &signal = signals[i];
//            parlay::sequence<double> calibrated_signal(signal.begin(), signal.end());
//            tstat_segmenter_t segmenter;
//            auto events = generate_events(calibrated_signal, segmenter);
//            auto quantized = quantize_signal(events);
//            parlay::sequence<u4> keys(quantized.size() - k + 1);
//            for (int j = 0; j < quantized.size() - k + 1; ++j)
//                keys[i] = encode_kmer(quantized.data() + j, k, 16, encode_qsig);
//            const auto [trg_id, trg_pos, support] = index.search(keys.data(), keys.size());
//            if (trg_id == -1) {
//                alignments[i].trg_name = '*', alignments[i].trg_pos = -1, alignments[i].presence = 0;
//            } else {
//                alignments[i].trg_name = index.headers[trg_id], alignments[i].trg_pos = trg_pos, alignments[i].presence = support;
//            }
//        });
//        return alignments;
//    }
//
//    std::vector<alignment_t> query_batch(std::vector<std::vector<int>> &signals, std::vector<double> &digitisations,
//                                         std::vector<double> &offsets, std::vector<double> &ranges) {
//        const auto B = signals.size();
//        std::vector<alignment_t> alignments(B);
//        parlay::for_each(parlay::iota(B), [&](size_t i){
//            auto &signal = signals[i];
//            parlay::sequence<double> calibrated_signal(signal.size());
//            for (int j = 0; j < signal.size(); ++j) {
//                calibrated_signal[j] = TO_PICOAMPS(signal[j], digitisations[i], offsets[i], ranges[i]);
//            }
//            tstat_segmenter_t segmenter;
//            auto events = generate_events(calibrated_signal, segmenter);
//            auto quantized = quantize_signal(events);
//            parlay::sequence<u4> keys(quantized.size() - k + 1);
//            for (int j = 0; j < quantized.size() - k + 1; ++j)
//                keys[i] = encode_kmer(quantized.data() + j, k, 16, encode_qsig);
//            const auto [trg_id, trg_pos, support] = index.search(keys.data(), keys.size());
//            if (trg_id == -1) {
//                alignments[i].trg_name = '*', alignments[i].trg_pos = -1, alignments[i].presence = 0;
//            } else {
//                alignments[i].trg_name = index.headers[trg_id], alignments[i].trg_pos = trg_pos, alignments[i].presence = support;
//            }
//        });
//        return alignments;
//    }
//
//
//};


struct Index {
    explicit Index(const std::string &input, int k=15, bool fwd_reverse=true, bool jaccard=true) {
        if (jaccard) idx = new j_index_t(k, 4, fwd_reverse);
        else idx = new c_index_t(k, 4, fwd_reverse);
        if (str_endswith(input.c_str(), ".cidx")) {
            // This is an index. Load it
            log_info("Loading index from %s", input.c_str());
            idx->load(input);
        } else if (str_endswith(input.c_str(), ".fa") || str_endswith(input.c_str(), ".fasta")) {
            log_info("Building index from %s", input.c_str());
            index_fasta(input.c_str(), idx);
        } else if (str_endswith(input.c_str(), ".fa.gz") || str_endswith(input.c_str(), ".fasta.gz")) {
            // zipped reference. will support it later
            log_error("This is not implemented yet");
        } else {
            log_error("Unknown input file format");
        }
        idx->init_query_buffers();
    }

    void dump(const std::string &basename) {
        std::string filename = basename + ".cidx";
        log_info("Dumping index to %s", filename.c_str());
        idx->dump(filename);
        log_info("Done.");
    }

    Alignment query(std::string &sequence) {
        auto result = idx->search(sequence);
        return {
            std::get<0>(result), std::get<1>(result),
                    static_cast<int>(std::get<2>(result)), std::get<3>(result), static_cast<int>(sequence.size())};
    }

    std::vector<Alignment> query_batch(const py::list& sequences) {
        auto nr = sequences.size();
        std::vector<Alignment> results(nr);
        parlay::for_each(parlay::iota(nr), [&](size_t i){
            auto sequence = sequences[i].cast<std::string>();
            results[i] = query(sequence);
        });
        return results;
    }



private:
    index_t *idx = nullptr;
};

// Binding the function to the Python module
PYBIND11_MODULE(_core, m) {
    py::class_<Alignment>(m, "Alignment")
            .def(py::init<>())  // Default constructor
            .def(py::init<const char*, bool, int, float, int>(),  // Parameterized constructor
                 py::arg("header"), py::arg("fwd"), py::arg("start"),
                 py::arg("pres_frac"), py::arg("qry_len"))
            .def_readonly("ctg", &Alignment::ctg)
            .def_readonly("r_st", &Alignment::r_st)
            .def_readonly("r_en", &Alignment::r_en)
            .def_readonly("strand", &Alignment::strand)
            .def_readonly("pres_frac", &Alignment::pres_frac);

    py::class_<Index>(m, "Index")
            .def(py::init<const std::string &, int, bool, bool>(),
                 py::arg("input"), py::arg("k")=15, py::arg("fwd_reverse")=true, py::arg("jaccard")=true)
            .def("dump", &Index::dump)
            .def("query", &Index::query)
            .def("query_batch", &Index::query_batch);
}