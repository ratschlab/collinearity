//
// Created by Sayan Goswami on 16.01.2025.
//

#include "collinearity.h"
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "argparse/argparse.hpp"

namespace py = pybind11;

static vector<const char*> kwargs_to_argv(const py::args& args, const py::kwargs& kwargs,
                                               vector<string>& arg_strings)  // Store actual argument strings
{
    vector<const char*> argv;  // Store pointers to C-style strings

    // First argument is conventionally the program name
    arg_strings.emplace_back("program_name");
    argv.push_back(arg_strings.back().c_str());

    // Handle positional arguments (flags)
    for (auto& arg : args) {
        string flag = "--" + arg.cast<string>();  // Convert to "--flag" format
        arg_strings.emplace_back(flag);
    }

    // Handle keyword arguments ("--key=value")
    for (auto& item : kwargs) {
        auto key = item.first.cast<string>();
        string value = py::str(item.second);  // Convert value to string
        if (key.size()==1) arg_strings.emplace_back("-" + key + "=" + value);
        else arg_strings.emplace_back("--" + key + "=" + value);
    }

    // Print for demonstration
//    for (auto &s : arg_strings)
//        cout << s << '\n';
//    printf("-----------\n");

    for (auto &s: arg_strings)
        argv.push_back(s.c_str());

    // Now `argc` and `argv.data()` can be used in a function expecting C-style args
    return argv;
}

struct rf_config_t : config_t {
    static rf_config_t init(int argc, char *argv[]) {
        auto c = argparse::parse<rf_config_t>(argc, argv);
        c.print();
        return c;
    }
};

struct Alignment {
    string ctg;
    int r_st = 0, r_en = 0, strand = 1;
    float pres_frac = 0.0f;

    Alignment() = default;

    Alignment(const char *header, bool fwd, int start, float pres_frac, int qry_len):
        ctg(header), r_st(start), r_en(start + qry_len), strand(fwd?1:-1), pres_frac(pres_frac) {}
};

struct Request {
    int channel = 0;
    string id;
    string seq;
    Request() = default;
    Request(int channel, string &id, string &seq): channel(channel), id(id), seq(seq) {}
};

struct Response {
    int channel = 0;
    string id;
    Alignment alignment;
    Response() = default;
    Response(int channel, string &id, Alignment &alignment): channel(channel), id(id), alignment(alignment) {}
};

struct ResponseGenerator {
    parlay::sequence<Response> responses;
    explicit ResponseGenerator(parlay::sequence<Response> &responses): responses(responses) {}
    Response next() {
        if (responses.empty())
            throw py::stop_iteration();
        else {
            auto response = responses.back();
            responses.pop_back();
            return response;
        }
    }
    ResponseGenerator& iter() {
        return *this;
    }
};

struct Index {
    explicit Index(string &input, const py::args& args, const py::kwargs& kwargs) {
        auto argv = kwargs_to_argv(args, kwargs, argvec);
        auto config = argparse::parse<rf_config_t>(argv.size(), argv.data());
        config.print();
        if (config._sort_block_size.empty()) config.sort_block_size = MEMPOOL_BLOCKSZ;
        else config.sort_block_size = hmsize2bytes(config._sort_block_size);
        config.sort_block_size = (config.sort_block_size < MEMPOOL_BLOCKSZ)? MEMPOOL_BLOCKSZ : config.sort_block_size;

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
            // gzipped reference. will support it later - todo
            log_error("This is not implemented yet");
        } else {
            log_error("Unknown input file format");
        }
        idx->init_query_buffers();
    }

    void dump(const string &basename) {
        string filename = basename + ".cidx";
        log_info("Dumping index to %s", filename.c_str());
        idx->dump(filename);
        log_info("Done.");
    }

    void load(const string &basename) {
        string filename = basename + ".cidx";
        log_info("Loading index from %s", filename.c_str());
        idx->load(filename);
        log_info("Done.");
    }

    Alignment query(string &sequence) {
        auto result = idx->search(sequence);
        return {
            get<0>(result), get<1>(result),
                    static_cast<int>(get<2>(result)), get<3>(result), static_cast<int>(sequence.size())};
    }

    vector<Alignment> query_batch(const py::list& sequences) {
        auto nr = sequences.size();
        vector<Alignment> results(nr);
        parlay::for_each(parlay::iota(nr), [&](size_t i){
            auto sequence = sequences[i].cast<string>();
            results[i] = query(sequence);
        });
        return results;
    }

    ResponseGenerator query_stream(const py::iterator& reads) {
        if (!stream_ready) log_error("Query stream is not ready.");
        parlay::sequence<Request> requests;
        for (auto &read: reads) {
            auto request = read.cast<Request>();
            requests.push_back(request);
        }
        auto responses = parlay::tabulate(requests.size(), [&](size_t i) {
            auto alignment = query(requests[i].seq);
            return Response(requests[i].channel, requests[i].id, alignment);
        });
        return ResponseGenerator(responses);
    }

protected:
    index_t *idx = nullptr;
    vector<string> argvec;
    bool stream_ready = false;
};

struct DynIndex {
    DynIndex(const py::args& args, const py::kwargs& kwargs) {
        auto argv = kwargs_to_argv(args, kwargs, argvec);
        auto config = argparse::parse<config_t>(argv.size(), argv.data());
        idx = new dindex_t(config);
    }

    void add(string &name, string &seq) {
        idx->add(name, parlay::make_slice(seq.data(), seq.data() + seq.size()));
    }

    void add_batch(const py::list& names, const py::list& sequences) { // maybe I'll find a better way later
        auto nr = sequences.size();
        for (int i = 0; i < nr; ++i) {
            auto name = names[i].cast<string>();
            auto seq = sequences[i].cast<string>();
            idx->add(name, parlay::make_slice(seq.data(), seq.data() + seq.size()));
        }
    }

    void merge() { idx->merge(); }

    Alignment query(string &sequence) {
        auto result = idx->search(sequence);
        return {
                get<0>(result), get<1>(result),
                static_cast<int>(get<2>(result)), get<3>(result), static_cast<int>(sequence.size())};
    }

    vector<Alignment> query_batch(const py::list& sequences) {
        auto nr = sequences.size();
        vector<Alignment> results(nr);
        parlay::for_each(parlay::iota(nr), [&](size_t i){
            auto sequence = sequences[i].cast<string>();
            results[i] = query(sequence);
        });
        return results;
    }

private:
    dindex_t *idx = nullptr;
    vector<string> argvec;
};


//class StreamProcessor {
//public:
//    explicit StreamProcessor(int n_threads = 1) {
//        for (unsigned int i = 0; i < n_threads; ++i) {
//            workers.emplace_back([this, i]() { this->worker_loop(i); });
//        }
//        ready = true;
//    }
//    ~StreamProcessor() {
//        if (ready) teardown();
//    }
//
//    void teardown() {
//        if (ready) {
//            const Request stop(-1, (string &) "", (string &) "");
//            for (auto &t : workers)
//                requests.enqueue(stop);
//            for (auto &t : workers)
//                if (t.joinable()) t.join();
//        }
//        ready = false;
//    }
//
//    ResponseGenerator batch_align(const py::iterator& reads) {
//        for (auto &read: reads) {
//            auto request = read.cast<Request>();
//            requests.enqueue(request);
//        }
//        return ResponseGenerator(responses);
//    }
//
//protected:
//    bool ready;
//    vector<thread> workers;
//    mc::ConcurrentQueue<Request> requests;
//    mc::ConcurrentQueue<Response> responses;
//
//    void worker_loop(u4 i) {
//        const u4 id = i;
//        bool done = false;
//        while (!done) {
//            Request request;
//            if (requests.try_dequeue(request)) {
//                if (request.channel < 0) done = true;
//                else {
//                    auto seq = request.seq;
//                    // process request - todo
//                    auto alignment = Alignment("", true, 0, 0.0f, 0);
//                    auto response = Response(request.channel, request.id, alignment);
//                    responses.enqueue(response);
//                    sitrep("%u - %s", id, seq.substr(0, 10).c_str());
//                }
//            } else this_thread::sleep_for(chrono::milliseconds(10));
//        }
//    }
//};


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
            .def(py::init<string&, const py::args&, const py::kwargs&>(), py::arg("input"))
            .def("dump", &Index::dump)
            .def("load", &Index::load)
            .def("query", &Index::query)
            .def("query_batch", &Index::query_batch)
            .def("query_stream", &Index::query_stream);

    py::class_<DynIndex>(m, "DynIndex")
            .def(py::init<const py::args&, const py::kwargs&>())
            .def("add", &DynIndex::add)
            .def("add_batch", &DynIndex::add_batch)
            .def("merge", &DynIndex::merge)
            .def("query", &DynIndex::query)
            .def("query_batch", &DynIndex::query_batch);

    py::class_<Request>(m, "Request")
            .def(py::init<int, string&, string&>(), py::arg("channel"), py::arg("id"), py::arg("seq"))
            .def_readwrite("channel", &Request::channel)
            .def_readwrite("id", &Request::id)
            .def_readwrite("seq", &Request::seq);

    py::class_<Response>(m, "Response")
            .def(py::init<int, string&, Alignment&>(), py::arg("channel"), py::arg("id"), py::arg("alignment"))
            .def_readwrite("channel", &Response::channel)
            .def_readwrite("id", &Response::id)
            .def_readwrite("alignment", &Response::alignment);

    py::class_<ResponseGenerator>(m, "ResponseGenerator")
            .def("__iter__", &ResponseGenerator::iter)
            .def("__next__", &ResponseGenerator::next);
}