//
// Created by Sayan Goswami on 09.12.2024.
//

#include "collinearity.h"
#include "tsl/hopscotch_map.h"

#define SANITY_CHECKS 1

using namespace klibpp;

static void align(std::vector<std::string> &qry_headers, std::vector<u4> &lengths, parlay::sequence<u4> &qkeys,
                  std::vector<std::string> &trg_headers, index_t *index,
                  tsl::hopscotch_map<u8, u4> *intercept_counts);
static inline void print_alignments(std::vector<std::string> &qry_headers, std::vector<std::string> &trg_headers,
                                    parlay::sequence<u4> &trg_ids, parlay::sequence<u8> &trg_posns);

void query(const char *filename, int k, int sigma, const size_t batch_sz, index_t *index, std::vector<std::string> &refnames) {
    // 1. read sequences in a batch
    // 2. create key-value pairs
    // 3. find collinear chains

//    idx_slice_t isl;
    tsl::hopscotch_map<u8, u4> intercept_counts[batch_sz];

    std::vector<std::string> headers;
    std::vector<u4> lengths;
    headers.reserve(batch_sz);
    lengths.reserve(batch_sz);
    parlay::sequence<u4> qKeys;
    // read sequences and convert to key-value pairs
    KSeq record;
    auto fd = open(filename, O_RDONLY);
    if (fd < 0) error("Could not open %s because %s.", filename, strerror(errno));
    auto ks = make_kstream(fd, read, mode::in);
    info("Begin query..");

    int nr = 0;
    while (ks >> record) {
        if (record.seq.size() > k) {
            headers.push_back(record.name);
            auto kmers = create_kmers(record.seq, k, sigma);
            qKeys.append(kmers.begin(), kmers.end());
            lengths.push_back(kmers.size());
            nr++;
            if (nr == batch_sz) {
                // query
                expect(headers.size() == nr);
                expect(lengths.size() == nr);
                align(headers, lengths, qKeys, refnames, index, intercept_counts);
                info("%d", nr);
                nr = 0;
                qKeys.clear();
                headers.clear();
                lengths.clear();
            }
        }
    }
    if (nr) {
        // query
        expect(headers.size() == nr);
        expect(lengths.size() == nr);
        align(headers, lengths, qKeys, refnames, index, intercept_counts);
        info("%d", nr);
        nr = 0;
        qKeys.clear();
        headers.clear();
        lengths.clear();
    }
    printf("\n");
    close(fd);
    info("Done.");
}

static void align(std::vector<std::string> &qry_headers, std::vector<u4> &lengths, parlay::sequence<u4> &qkeys,
                  std::vector<std::string> &trg_headers, index_t *index,
                  tsl::hopscotch_map<u8, u4> *intercept_counts) {
    const auto B = qry_headers.size();
    expect(B == lengths.size());
    auto qry_offsets = parlay::scan(lengths);
    parlay::sequence<u4> trg_ref_id(B);
    parlay::sequence<u8> trg_pos(B);

//    parlay::for_each(parlay::iota(B), [&](size_t i){
    for (u8 i = 0; i < B; ++i){
        intercept_counts[i].clear();
        auto qry_offset = qry_offsets.first[i];
        auto qry_size = lengths[i];
        for (u4 qry_pos = 0, j = qry_offset; j < qry_offset + qry_size; ++j, ++qry_pos) {
            auto &values = index->get(qkeys[j]);
            if (!values.empty()) {
                for (auto v : values) {
                    u8 ref_id = get_id_from(v);
                    u8 ref_pos = get_pos_from(v);
                    uint intercept = (ref_pos > qry_pos) ? (ref_pos - qry_pos) : 0;
                    intercept /= bandwidth;
                    u8 key = make_key_from(ref_id, intercept);
                    if (intercept_counts[i].contains(key)) intercept_counts[i][key]++;
                    else intercept_counts[i][key] = 1;
                    if (intercept >= bandwidth) {
                        intercept -= bandwidth;
                        key = make_key_from(ref_id, intercept);
                        if (intercept_counts[i].contains(key)) intercept_counts[i][key]++;
                        else intercept_counts[i][key] = 1;
                    }
                }
            }
//            auto [start, end] = isl.get(qkeys[j]);
//            for (auto p = start; p < end; ++p) {
//                auto v = *p;
//                if (v) {
//                    u8 ref_id = get_id_from(v);
//                    u8 ref_pos = get_pos_from(v);
//                    uint intercept = (ref_pos > qry_pos) ? (ref_pos - qry_pos) : 0;
//                    intercept /= bandwidth;
//                    u8 key = make_key_from(ref_id, intercept);
//                    if (intercept_counts[i].contains(key)) intercept_counts[i][key]++;
//                    else intercept_counts[i][key] = 1;
//                    if (intercept >= bandwidth) {
//                        intercept -= bandwidth;
//                        key = make_key_from(ref_id, intercept);
//                        if (intercept_counts[i].contains(key)) intercept_counts[i][key]++;
//                        else intercept_counts[i][key] = 1;
//                    }
//                }
//            }
        }
        // traverse through the map and check the intercept with the maximum votes
        uint max_nvotes = 0;
        uint64_t best_key = 0;
        for(const auto& kv : intercept_counts[i]) {
            if (kv.second > max_nvotes) {
                max_nvotes = kv.second, best_key = kv.first;
            }
        }
        // in this case, discovery_fraction is the frac. of kmers that support an intercept
        if ((max_nvotes * 1.0) / qry_size >= presence_fraction) {
            trg_ref_id[i] = get_id_from(best_key);
            trg_pos[i] = get_pos_from(best_key) * bandwidth;
        } else {
            trg_ref_id[i] = -1;
            trg_pos[i] = -1;
        }
    }
//    );
    print_alignments(qry_headers, trg_headers, trg_ref_id, trg_pos);
}

static inline void print_alignments(std::vector<std::string> &qry_headers, std::vector<std::string> &trg_headers,
                             parlay::sequence<u4> &trg_ids, parlay::sequence<u8> &trg_posns) {
    for (u4 i = 0; i < qry_headers.size(); ++i){
        if (trg_ids[i] == -1 && trg_posns[i] == -1)
            printf("%s\t*\t0\n", qry_headers[i].c_str());
        else
            printf("%s\t%s\t%ld\n", qry_headers[i].c_str(), trg_headers[trg_ids[i]].c_str(), trg_posns[i]);
    }
}