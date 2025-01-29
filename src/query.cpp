//
// Created by Sayan Goswami on 09.12.2024.
//

#include "collinearity.h"
#include "../external/slow5lib/include/slow5/slow5.h"

#define SANITY_CHECKS 1

using namespace klibpp;

static void align(const std::vector<std::string> &qry_headers, const std::vector<u4> &lengths,
                  const parlay::sequence<u4> &qkeys, index_t &index);
static void align_raw(const std::vector<std::string> &qry_headers, const std::vector<u4> &lengths,
                      const std::vector<double> &digitizations, const std::vector<double> &offsets,
                      const std::vector<double> &ranges,
                      parlay::sequence<int16_t > &raw_signals,
                      index_t &index, int k, int sigma);
static inline void print_alignments(const std::vector<std::string> &qry_headers, const std::vector<std::string> &trg_headers,
                                    const parlay::sequence<u4> &trg_ids, const parlay::sequence<u8> &trg_posns);

void query(const char *filename, int k, int sigma, const size_t batch_sz, index_t &index) {
    // 1. read sequences in a batch
    // 2. create key-value pairs
    // 3. find collinear chains

    index.init_query_buffers();
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
            auto kmers = create_kmers(record.seq, k, sigma, encode_dna);
            qKeys.append(kmers.begin(), kmers.end());
            lengths.push_back(kmers.size());
            nr++;
            if (nr == batch_sz) {
                // query
                expect(headers.size() == nr);
                expect(lengths.size() == nr);
                align(headers, lengths, qKeys, index);
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
        align(headers, lengths, qKeys, index);
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

void query_raw(const char *filename, int k, int sigma, const size_t batch_sz, index_t &index) {
    index.init_query_buffers();
    std::vector<std::string> headers;
    std::vector<u4> lengths;
    std::vector<double> digitizations, offsets, ranges;
    headers.reserve(batch_sz);
    lengths.reserve(batch_sz);
    digitizations.reserve(batch_sz), offsets.reserve(batch_sz), ranges.reserve(batch_sz);

    parlay::sequence<int16_t> raw_signals;

    int nr = 0;
    u8 total_nr = 0;

    auto sp = slow5_open(filename, "r");
    if (!sp) error("Could not open file %s", filename);

    slow5_rec_t *rec = nullptr; //slow5 record to be read
    int ret=0; //for return value

    //iterate through the file until end
    while((ret = slow5_get_next(&rec,sp)) >= 0){
        headers.emplace_back(rec->read_id);
        lengths.push_back(rec->len_raw_signal);
        digitizations.push_back(rec->digitisation);
        offsets.push_back(rec->offset);
        ranges.push_back(rec->range);
        raw_signals.append(rec->raw_signal, rec->raw_signal + rec->len_raw_signal);
        //double pA = TO_PICOAMPS(rec->raw_signal[i],rec->digitisation,rec->offset,rec->range);
        nr++;
        if (nr == batch_sz) {
            // query
            expect(headers.size() == nr);
            expect(lengths.size() == nr);
            align_raw(headers, lengths, digitizations, offsets, ranges, raw_signals, index, k, sigma);
            total_nr += nr;
            sitrep("%lu", total_nr);
            nr = 0;
            headers.clear(), lengths.clear(); digitizations.clear(), ranges.clear(), offsets.clear(), raw_signals.clear();
        }
    }

    if (nr) {
        // query
        expect(headers.size() == nr);
        expect(lengths.size() == nr);
        align_raw(headers, lengths, digitizations, offsets, ranges, raw_signals, index, k, sigma);
        total_nr += nr;
        sitrep("%lu", total_nr);
        nr = 0;
        headers.clear(), lengths.clear(); digitizations.clear(), ranges.clear(), offsets.clear(), raw_signals.clear();
    }
    stderrflush;

    if(ret != SLOW5_ERR_EOF)  //check if proper end of file has been reached
        error("Error in slow5_get_next. Error code %d\n",ret);

    //free the SLOW5 record
    slow5_rec_free(rec);

    //close the SLOW5 file
    slow5_close(sp);
}

static void align(const std::vector<std::string> &qry_headers, const std::vector<u4> &lengths,
                  const parlay::sequence<u4> &qkeys, index_t &index) {
    const auto B = qry_headers.size();
    expect(B == lengths.size());
    auto qry_offsets = parlay::scan(lengths);
    parlay::sequence<u4> trg_ref_id(B);
    parlay::sequence<u8> trg_pos(B);
    parlay::sequence<float> presence(B);

    parlay::for_each(parlay::iota(B), [&](size_t i){
        auto qry_offset = qry_offsets.first[i];
        auto qry_size = lengths[i];
        std::tie(trg_ref_id[i], trg_pos[i], presence[i]) = index.search(qkeys.data() + qry_offset, qry_size);
    });
    print_alignments(qry_headers, index.headers, trg_ref_id, trg_pos);
}

static inline void print_alignments(const std::vector<std::string> &qry_headers, const std::vector<std::string> &trg_headers,
                             const parlay::sequence<u4> &trg_ids, const parlay::sequence<u8> &trg_posns) {
    for (u4 i = 0; i < qry_headers.size(); ++i){
        if (trg_ids[i] == (u4)-1)
            printf("%s\t*\t0\n", qry_headers[i].c_str());
        else
            printf("%s\t%s\t%ld\n", qry_headers[i].c_str(), trg_headers[trg_ids[i]].c_str(), trg_posns[i]);
    }
}

static void align_raw(const std::vector<std::string> &qry_headers, const std::vector<u4> &lengths,
                      const std::vector<double> &digitizations, const std::vector<double> &offsets,
                      const std::vector<double> &ranges,
                      parlay::sequence<int16_t > &raw_signals,
                      index_t &index, const int k, const int sigma) {
    const auto B = qry_headers.size();
    expect(B == lengths.size());
    const auto [qry_offsets, total_qry_len] = parlay::scan(lengths);
    parlay::sequence<u4> trg_ref_id(B);
    parlay::sequence<u8> trg_pos(B);
    parlay::sequence<float> presence(B);

    parlay::for_each(parlay::iota(B), [&](size_t i) {
        auto qry_offset = qry_offsets[i];
        auto signal_length = lengths[i];
        auto signal_digitization = digitizations[i];
        auto signal_offset = offsets[i];
        auto signal_range = ranges[i];
        parlay::sequence<double> calibrated_signal(signal_length);
        for (int j = 0; j < signal_length; ++j) {
            calibrated_signal[j] = TO_PICOAMPS(raw_signals[qry_offset + j], signal_digitization, signal_offset,
                                               signal_range);
        }
        tstat_segmenter_t segmenter;
        auto events = generate_events(calibrated_signal, segmenter);
        auto quantized = quantize_signal(events);
        parlay::sequence<u4> keys(quantized.size() - k + 1);
        for (int j = 0; j < quantized.size() - k + 1; ++j)
            keys[j] = encode_kmer(quantized.data() + j, k, sigma, encode_qsig);
        std::tie(trg_ref_id[i], trg_pos[i], presence[i]) = index.search(keys.data(), keys.size());
    });
    print_alignments(qry_headers, index.headers, trg_ref_id, trg_pos);
}