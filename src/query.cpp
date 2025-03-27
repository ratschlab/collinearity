//
// Created by Sayan Goswami on 09.12.2024.
//

#include "collinearity.h"
#include "../external/slow5lib/include/slow5/slow5.h"

#define SANITY_CHECKS 1

using namespace klibpp;

const char STRAND[2] = {'-', '+'};

std::vector<std::tuple<std::string, bool, u4, float>> query_batch(index_t *idx, std::vector<std::string> &sequences) {
    auto nr = sequences.size();
    std::vector<std::tuple<std::string, bool, u4, float>> results;
    parlay::for_each(parlay::iota(nr), [&](size_t i){
        auto result = idx->search(sequences[i]);
        return std::make_tuple(
                std::string(std::get<0>(result)), std::get<1>(result), std::get<2>(result), std::get<3>(result));
    });
    return results;
}

void query_fasta(index_t *idx, const char* fasta_filename, int batch_sz) {
    idx->init_query_buffers();
    std::vector<std::string> headers, sequences;
    headers.reserve(batch_sz);
    sequences.reserve(batch_sz);

    // read sequences and convert to key-value pairs
    KSeq record;
    auto fd = open(fasta_filename, O_RDONLY);
    if (fd < 0) log_error("Could not open %s because %s.", fasta_filename, strerror(errno));
    auto ks = make_kstream(fd, read, mode::in);
    log_info("Begin query..");
    int nr = 0;
    u8 total_nr = 0;
    while (ks >> record) {
        headers.emplace_back(record.name);
        sequences.emplace_back(record.seq);
        nr++;
        if (nr == batch_sz) {
            auto results = parlay::tabulate(nr, [&](size_t i) {
                return idx->search(sequences[i]);
            });
            for (u4 i = 0; i < nr; ++i)
                printf("%s\t%s\t%c\t%d\t%f\n", headers[i].c_str(), std::get<0>(results[i]),
                       STRAND[(int)std::get<1>(results[i])], std::get<2>(results[i]), std::get<3>(results[i]));
            total_nr += nr;
            sitrep("%lu", total_nr);
            nr = 0;
            headers.clear();
            sequences.clear();
        }
    }
    if (nr) {
        auto results = parlay::tabulate(nr, [&](size_t i) {
            return idx->search(sequences[i]);
        });
        for (u4 i = 0; i < nr; ++i)
            printf("%s\t%s\t%c\t%d\t%f\n", headers[i].c_str(), std::get<0>(results[i]),
                   STRAND[(int)std::get<1>(results[i])], std::get<2>(results[i]), std::get<3>(results[i]));
        total_nr += nr;
        sitrep("%lu", total_nr);
        nr = 0;
        headers.clear();
        sequences.clear();
    }
    stderrflush;
    close(fd);
    log_info("Done.");
}

void query_blow5(index_t *idx, const char* blow5_filename, int batch_sz) {
    throw "Not implemented";
}

//static void align(const std::vector<std::string> &qry_headers, const std::vector<u4> &lengths,
//                  const parlay::sequence<u4> &qkeys, sindex_t &index);
//static void align_raw(const std::vector<std::string> &qry_headers, const std::vector<u4> &lengths,
//                      const std::vector<double> &digitizations, const std::vector<double> &offsets,
//                      const std::vector<double> &ranges,
//                      parlay::sequence<int16_t > &raw_signals,
//                      sindex_t &index, int k, int sigma);
//static inline void print_alignments(const std::vector<std::string> &qry_headers, const std::vector<std::string> &trg_headers,
//                                    const parlay::sequence<u4> &trg_ids, const parlay::sequence<u8> &trg_posns);
//
//void query(const char *filename, int k, int sigma, const size_t batch_sz, sindex_t &index) {
//    // 1. read sequences in a batch
//    // 2. create key-value pairs
//    // 3. find collinear chains
//
//    index.init_query_buffers();
//    std::vector<std::string> headers;
//    std::vector<u4> lengths;
//    headers.reserve(batch_sz);
//    lengths.reserve(batch_sz);
//    parlay::sequence<u4> qKeys;
//    // read sequences and convert to key-value pairs
//    KSeq record;
//    auto fd = open(filename, O_RDONLY);
//    if (fd < 0) log_error("Could not open %s because %s.", filename, strerror(errno));
//    auto ks = make_kstream(fd, read, mode::in);
//    log_info("Begin query..");
//    int nr = 0;
//    u8 total_nr = 0;
//    while (ks >> record) {
//        if (record.seq.size() > k) {
//            headers.push_back(record.name);
//            auto kmers = create_kmers(record.seq, k, sigma, encode_dna);
//            qKeys.append(kmers.begin(), kmers.end());
//            lengths.push_back(kmers.size());
//            nr++;
//            if (nr == batch_sz) {
//                // query
//                expect(headers.size() == nr);
//                expect(lengths.size() == nr);
//                align(headers, lengths, qKeys, index);
//                total_nr += nr;
//                sitrep("%lu", total_nr);
//                nr = 0;
//                qKeys.clear();
//                headers.clear();
//                lengths.clear();
//            }
//        }
//    }
//    if (nr) {
//        // query
//        expect(headers.size() == nr);
//        expect(lengths.size() == nr);
//        align(headers, lengths, qKeys, index);
//        total_nr += nr;
//        sitrep("%lu", total_nr);
//        nr = 0;
//        qKeys.clear();
//        headers.clear();
//        lengths.clear();
//    }
//    stderrflush;
//    close(fd);
//    log_info("Done.");
//}
//
//void query_jaccard(const char *filename, int k, int sigma, size_t batch_sz, jindex_t &index) {
//    index.init_query_buffers();
//    std::vector<std::string> headers, sequences;
//    headers.reserve(batch_sz);
//    sequences.reserve(batch_sz);
//    // read sequences and convert to key-value pairs
//    KSeq record;
//    auto fd = open(filename, O_RDONLY);
//    if (fd < 0) log_error("Could not open %s because %s.", filename, strerror(errno));
//    auto ks = make_kstream(fd, read, mode::in);
//    log_info("Begin query..");
//    int nr = 0;
//    u8 total_nr = 0;
//    while (ks >> record) {
//        if (record.seq.size() > k) {
//            headers.push_back(record.name);
//            sequences.emplace_back(record.seq);
//            nr++;
//            if (nr == batch_sz) {
//                auto results = parlay::tabulate(nr, [&](size_t i) {
//                    return index.get(sequences[i]);
//                });
//                for (u4 i = 0; i < nr; ++i)
//                    printf("%s\t%s\t%d\n", headers[i].c_str(), get<0>(results[i]), get<1>(results[i]));
//                total_nr += nr;
//                sitrep("%lu", total_nr);
//                nr = 0;
//                headers.clear();
//                sequences.clear();
//            }
//        }
//    }
//    if (nr) {
//        auto results = parlay::tabulate(nr, [&](size_t i) {
//            return index.get(sequences[i]);
//        });
//        for (u4 i = 0; i < nr; ++i)
//            printf("%s\t%s\t%d\n", headers[i].c_str(), get<0>(results[i]), get<1>(results[i]));
//        total_nr += nr;
//        sitrep("%lu", total_nr);
//        nr = 0;
//        headers.clear();
//        sequences.clear();
//    }
//}
//
//void query_raw(const char *filename, int k, int sigma, const size_t batch_sz, sindex_t &index) {
//    index.init_query_buffers();
//    std::vector<std::string> headers;
//    std::vector<u4> lengths;
//    std::vector<double> digitizations, offsets, ranges;
//    headers.reserve(batch_sz);
//    lengths.reserve(batch_sz);
//    digitizations.reserve(batch_sz), offsets.reserve(batch_sz), ranges.reserve(batch_sz);
//
//    parlay::sequence<int16_t> raw_signals;
//
//    int nr = 0;
//    u8 total_nr = 0;
//
//    auto sp = slow5_open(filename, "r");
//    if (!sp) log_error("Could not open file %s", filename);
//
//    slow5_rec_t *rec = nullptr; //slow5 record to be read
//    int ret=0; //for return value
//
//    log_info("Querying..");
//
//    //iterate through the file until end
//    while((ret = slow5_get_next(&rec,sp)) >= 0){
//        headers.emplace_back(rec->read_id);
//        lengths.push_back(rec->len_raw_signal);
//        digitizations.push_back(rec->digitisation);
//        offsets.push_back(rec->offset);
//        ranges.push_back(rec->range);
//        raw_signals.append(rec->raw_signal, rec->raw_signal + rec->len_raw_signal);
//        //double pA = TO_PICOAMPS(rec->raw_signal[i],rec->digitisation,rec->offset,rec->range);
//        nr++;
//        if (nr == batch_sz) {
//            sitrep("----------");
//            // query
//            expect(headers.size() == nr);
//            expect(lengths.size() == nr);
//            align_raw(headers, lengths, digitizations, offsets, ranges, raw_signals, index, k, sigma);
//            total_nr += nr;
//            sitrep("%lu", total_nr);
//            nr = 0;
//            headers.clear(), lengths.clear(); digitizations.clear(), ranges.clear(), offsets.clear(), raw_signals.clear();
//        }
//    }
//
//    if (nr) {
//        // query
//        expect(headers.size() == nr);
//        expect(lengths.size() == nr);
//        align_raw(headers, lengths, digitizations, offsets, ranges, raw_signals, index, k, sigma);
//        total_nr += nr;
//        sitrep("%lu", total_nr);
//        nr = 0;
//        headers.clear(), lengths.clear(); digitizations.clear(), ranges.clear(), offsets.clear(), raw_signals.clear();
//    }
//    stderrflush;
//
//    if(ret != SLOW5_ERR_EOF)  //check if proper end of file has been reached
//        log_error("Error in slow5_get_next. Error code %d\n",ret);
//
//    //free the SLOW5 record
//    slow5_rec_free(rec);
//
//    //close the SLOW5 file
//    slow5_close(sp);
//    log_info("Done.");
//}
//
//static void align(const std::vector<std::string> &qry_headers, const std::vector<u4> &lengths,
//                  const parlay::sequence<u4> &qkeys, sindex_t &index) {
//    const auto B = qry_headers.size();
//    expect(B == lengths.size());
//    auto qry_offsets = parlay::scan(lengths);
//    parlay::sequence<u4> trg_ref_id(B);
//    parlay::sequence<u8> trg_pos(B);
//    parlay::sequence<float> presence(B);
//
//    parlay::for_each(parlay::iota(B), [&](size_t i){
//        auto qry_offset = qry_offsets.first[i];
//        auto qry_size = lengths[i];
//        std::tie(trg_ref_id[i], trg_pos[i], presence[i]) = index.search(qkeys.data() + qry_offset, qry_size);
//    });
//    print_alignments(qry_headers, index.headers, trg_ref_id, trg_pos);
//}
//
//static inline void print_alignments(const std::vector<std::string> &qry_headers, const std::vector<std::string> &trg_headers,
//                             const parlay::sequence<u4> &trg_ids, const parlay::sequence<u8> &trg_posns) {
//    for (u4 i = 0; i < qry_headers.size(); ++i){
//        if (trg_ids[i] == (u4)-1)
//            printf("%s\t*\t0\n", qry_headers[i].c_str());
//        else
//            printf("%s\t%s\t%ld\n", qry_headers[i].c_str(), trg_headers[trg_ids[i]].c_str(), trg_posns[i]);
//    }
//}
//
//static void align_raw(const std::vector<std::string> &qry_headers, const std::vector<u4> &lengths,
//                      const std::vector<double> &digitizations, const std::vector<double> &offsets,
//                      const std::vector<double> &ranges,
//                      parlay::sequence<int16_t > &raw_signals,
//                      sindex_t &index, const int k, const int sigma) {
//    const auto B = qry_headers.size();
//    expect(B == lengths.size());
//    const auto [qry_offsets, total_qry_len] = parlay::scan(lengths);
//    parlay::sequence<u4> trg_ref_id(B);
//    parlay::sequence<u8> trg_pos(B);
//    parlay::sequence<float> presence(B);
//
//    parlay::for_each(parlay::iota(B), [&](size_t i) {
//        auto qry_offset = qry_offsets[i];
//        auto signal_length = lengths[i];
//        auto signal_digitization = digitizations[i];
//        auto signal_offset = offsets[i];
//        auto signal_range = ranges[i];
//        parlay::sequence<double> calibrated_signal(signal_length);
//        for (int j = 0; j < signal_length; ++j) {
//            calibrated_signal[j] = TO_PICOAMPS(raw_signals[qry_offset + j], signal_digitization, signal_offset,
//                                               signal_range);
//        }
//        tstat_segmenter_t segmenter;
//        auto events = generate_events(calibrated_signal, segmenter);
//        auto quantized = quantize_signal_simple(events);
//        parlay::sequence<u4> keys(quantized.size() - k + 1);
//        for (int j = 0; j < quantized.size() - k + 1; ++j)
//            keys[j] = encode_kmer(quantized.data() + j, k, sigma, encode_qsig);
//        std::tie(trg_ref_id[i], trg_pos[i], presence[i]) = index.search(keys.data(), keys.size());
//    });
//    print_alignments(qry_headers, index.headers, trg_ref_id, trg_pos);
//}