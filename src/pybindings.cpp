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
    index_t *index = nullptr;
    std::vector<std::string> refnames;
    heavyhitter_ht_t *hhs = nullptr;

    explicit rt_index_t(int batch_sz, const std::string &filename) {
        if (str_endswith(filename.c_str(), ".cidx"))
            std::tie(index, refnames) = load_index((char *)filename.c_str());
        else std::tie(index, refnames) = process_fasta(filename.c_str(), KMER_LENGTH, SIGMA);
        print_idx_info(index);
        hhs = new heavyhitter_ht_t[batch_sz];
    }

    void dump(const std::string &basename) {
        dump_index(index, refnames, basename.c_str());
    }

    std::vector<alignment_t> query_batch(std::vector<std::string> &sequences) {
        const auto B = sequences.size();
        std::vector<alignment_t> alignments(B);
        parlay::for_each(parlay::iota(B), [&](size_t i){
            auto &hh = hhs[i];
            const auto &sequence = sequences[i];
            const auto n = sequence.size();
            auto &alignment = alignments[i];
            hh.reset();
            u4 kmer = 0;
            for (int j = 0; j < KMER_LENGTH-1; ++j) kmer = kmer * SIGMA + encode_dna(sequence[j]);
            for (int j = 0; j < n - KMER_LENGTH + 1; j++) {
                kmer = kmer * SIGMA + encode_dna(sequence[KMER_LENGTH - 1 + j]);
                const auto &[vbegin, vend] = index->get(kmer);
                for (auto v = vbegin; v != vend; ++v) {
                    u8 ref_id = get_id_from(*v);
                    u8 ref_pos = get_pos_from(*v);
                    u8 intercept = (ref_pos > j) ? (ref_pos - j) : 0;
                    intercept /= bandwidth;
                    u8 key = make_key_from(ref_id, intercept);
                    hh.insert(key);
                    if (intercept >= bandwidth) {
                        intercept -= bandwidth;
                        key = make_key_from(ref_id, intercept);
                        hh.insert(key);
                    }
                }
            }

            // in this case, discovery_fraction is the frac. of kmers that support an intercept
            alignment.presence = (hh.top_count * 1.0f) / n;
            if (alignment.presence >= presence_fraction) {
                alignment.trg_name = refnames[get_id_from(hh.top_key)];
                alignment.trg_pos = get_pos_from(hh.top_key) * bandwidth;
            }
        });
        return alignments;
    }
};


// Binding the function to the Python module
PYBIND11_MODULE(_core, m) {

    pybind11::class_<rt_index_t>(m, "rt_index_t")
        .def(py::init<const int, const std::string &>())
        .def("dump", &rt_index_t::dump)
        .def("query_batch", &rt_index_t::query_batch);

    pybind11::class_<alignment_t>(m, "alignment_t")
            .def_readonly("trg_name", &alignment_t::trg_name)
            .def_readonly("trg_pos", &alignment_t::trg_pos)
            .def_readonly("presence", &alignment_t::presence);
}