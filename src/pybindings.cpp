//
// Created by Sayan Goswami on 16.01.2025.
//

#include "collinearity.h"
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "argparse/argparse.hpp"

namespace py = pybind11;

//static std::vector<const char*> kwargs_to_argv(const py::kwargs& kwargs, std::vector<std::string> &args) {
//    std::vector<const char*> argv;
//
//    // First argument is usually the program name, conventionally
//    args.emplace_back("program_name");
//    argv.push_back(args.back().c_str());
//
//    // Convert kwargs into "--key=value" format
//    for (auto& item : kwargs) {
//        std::string key = item.first.cast<std::string>();
//        std::string value = py::str(item.second);  // Convert value to string
//        args.emplace_back("--" + key + "=" + value);
//        argv.push_back(args.back().c_str());  // Store pointer to the new string
//    }
//
//    return argv;
//}

static std::vector<const char*> kwargs_to_argv(const py::args& args, const py::kwargs& kwargs,
                                               std::vector<std::string>& arg_strings)  // Store actual argument strings
{
    std::vector<const char*> argv;  // Store pointers to C-style strings

    // First argument is conventionally the program name
    arg_strings.emplace_back("program_name");
    argv.push_back(arg_strings.back().c_str());

    // Handle positional arguments (flags)
    for (auto& arg : args) {
        std::string flag = "--" + arg.cast<std::string>();  // Convert to "--flag" format
        arg_strings.emplace_back(flag);
    }

    // Handle keyword arguments ("--key=value")
    for (auto& item : kwargs) {
        std::string key = item.first.cast<std::string>();
        std::string value = py::str(item.second);  // Convert value to string
        if (key.size()==1) arg_strings.emplace_back("-" + key + "=" + value);
        else arg_strings.emplace_back("--" + key + "=" + value);
    }

    // Print for demonstration
//    for (auto &s : arg_strings)
//        std::cout << s << '\n';
//    printf("-----------\n");

    for (auto &s: arg_strings)
        argv.push_back(s.c_str());

    // Now `argc` and `argv.data()` can be used in a function expecting C-style args
    return argv;
}

struct Alignment {
    std::string ctg;
    int r_st, r_en, strand;
    float pres_frac;

    Alignment() = default;

    Alignment(const char *header, bool fwd, int start, float pres_frac, int qry_len):
        ctg(header), r_st(start), r_en(start + qry_len), strand(fwd?1:-1), pres_frac(pres_frac) {}
};

struct Index {
    explicit Index(std::string &input, const py::args& args, const py::kwargs& kwargs) {
        auto argv = kwargs_to_argv(args, kwargs, argvec);
        auto config = argparse::parse<config_t>(argv.size(), argv.data());
//        config.print();
        if (config.jaccard) log_info("Creating Jaccard index");
        else log_info("Creating a normal index");
        if (config.fwd_rev) log_info("Indexing both fwd and rev references");
        else log_info("Indexing fwd references only");

        if (config.jaccard) idx = new j_index_t(config);
        else idx = new c_index_t(config);

        if (str_endswith(input.c_str(), ".cidx")) {
            // This is an index. Load it
            log_info("Loading index from %s", input.c_str());
            idx->load(input);
        } else if (str_endswith(input.c_str(), ".fa") || str_endswith(input.c_str(), ".fasta")) {
            log_info("Building index from %s", input.c_str());
            index_fasta(input, idx);
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
    std::vector<std::string> argvec;
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
            .def(py::init<std::string&, const py::args&, const py::kwargs&>(), py::arg("input"))
            .def("dump", &Index::dump)
            .def("query", &Index::query)
            .def("query_batch", &Index::query_batch);
}