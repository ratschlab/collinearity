//
// Created by Sayan Goswami on 09.12.2024.
//

#include "collinearity.h"
#include "xxhash.h"
#include "hash_table8.hpp"

#define SANITY_CHECKS 1

using namespace klibpp;

// CMS estimate is within .8 percent of the count with a 99.5 percent probability
#define CMS_ERROR_RATE .008
#define CMS_CONFIDENCE .995
#define LOG_TWO 0.6931471805599453
// width = ceil(2 / error_rate)
#define CMS_WIDTH 1021
//const int CMS_WIDTH = (int)ceil(2.0/CMS_ERROR_RATE);
// depth = ceil((-1 * log(1 - confidence)) / LOG_TWO)
#define CMS_DEPTH 8
//const int CMS_DEPTH = (int)ceil((-1 * log(1 - CMS_CONFIDENCE)) / LOG_TWO);

struct [[maybe_unused]] heavyhitter_cms_t {
    u1 buf[CMS_DEPTH][CMS_WIDTH] = {0};
    u8 top_key = -1;
    u2 top_count = 0;
    void insert(const u8 key) {
        u2 count = -1;
        for (int i = 0; i < CMS_DEPTH; i+=2) {
            u8 hash = XXH64_hash64(key, i);
            u4 l1 = ((u4)hash) % CMS_WIDTH, l2 = ((u4)(hash>>32)) % CMS_WIDTH;
            buf[i][l1]++;
            buf[i+1][l2]++;
            count = MIN(count, buf[i][l1]);
            count = MIN(count, buf[i+1][l2]);
        }
        if (count > top_count)
            top_count = count, top_key = key;
    }
    void reset() { memset(buf, 0, sizeof(buf)), top_key=-1, top_count = 0; }
};

struct [[maybe_unused]] heavyhitter_ht_t {
    emhash8::HashMap<u8,u1> counts;
    u8 top_key = -1;
    u2 top_count = 0;
    void insert(const u8 key) {
        u4 count = counts[key]++;
        if (count > top_count)
            top_count = count, top_key = key;
    }
    void reset() { counts.clear(), top_key=-1, top_count = 0; }
};

static void align(std::vector<std::string> &qry_headers, std::vector<u4> &lengths, parlay::sequence<u4> &qkeys,
                  std::vector<std::string> &trg_headers, index_t *index, heavyhitter_ht_t *hhs);
static inline void print_alignments(std::vector<std::string> &qry_headers, std::vector<std::string> &trg_headers,
                                    parlay::sequence<u4> &trg_ids, parlay::sequence<u8> &trg_posns);

void query(const char *filename, int k, int sigma, const size_t batch_sz, index_t *index, std::vector<std::string> &refnames) {
    // 1. read sequences in a batch
    // 2. create key-value pairs
    // 3. find collinear chains

    heavyhitter_ht_t hhs[parlay::num_workers()];
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
    u8 total_nr = 0;
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
                align(headers, lengths, qKeys, refnames, index, hhs);
                total_nr += nr;
                sitrep("%lu", total_nr);
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
        align(headers, lengths, qKeys, refnames, index, hhs);
        total_nr += nr;
        sitrep("%lu", total_nr);
        nr = 0;
        qKeys.clear();
        headers.clear();
        lengths.clear();
    }
    stderrflush;
    close(fd);
    info("Done.");
}

static void align(std::vector<std::string> &qry_headers, std::vector<u4> &lengths, parlay::sequence<u4> &qkeys,
                  std::vector<std::string> &trg_headers, index_t *index, heavyhitter_ht_t *hhs) {
    const auto B = qry_headers.size();
    expect(B == lengths.size());
    auto qry_offsets = parlay::scan(lengths);
    parlay::sequence<u4> trg_ref_id(B);
    parlay::sequence<u8> trg_pos(B);

    parlay::for_each(parlay::iota(B), [&](size_t i){
        auto &hh = hhs[parlay::worker_id()];
        hh.reset();
        auto qry_offset = qry_offsets.first[i];
        auto qry_size = lengths[i];
        for (u4 qry_pos = 0, j = qry_offset; j < qry_offset + qry_size; ++j, ++qry_pos) {
            const auto &[vbegin, vend] = index->get(qkeys[j]);
            for (auto v = vbegin; v != vend; ++v) {
                u8 ref_id = get_id_from(*v);
                u8 ref_pos = get_pos_from(*v);
                u8 intercept = (ref_pos > qry_pos) ? (ref_pos - qry_pos) : 0;
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
        if ((hh.top_count * 1.0) / qry_size >= presence_fraction) {
            trg_ref_id[i] = get_id_from(hh.top_key);
            trg_pos[i] = get_pos_from(hh.top_key) * bandwidth;
        } else {
            trg_ref_id[i] = -1;
            trg_pos[i] = -1;
        }
    });
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