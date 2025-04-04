//
// Created by Sayan Goswami on 18.02.2025.
//
#include "collinearity.h"
#include "index.h"
#include "vbyte.h"


// I am assuming that the number of sequences won't exceed 2^20 (1M)
// and the longest sequence won't be longer than 2^40 (1T)
// if that's not the case, change the following line accordingly
#define ref_id_nbits 20
#define ref_len_nbits (64 - ref_id_nbits)
#define ref_id_bitmask (((1ULL) << ref_len_nbits) - 1)
#define make_key_from(id, pos) (((u8)(id)) << (ref_len_nbits) | (pos))
#define get_id_from(key) ((key) >> ref_len_nbits)
#define get_pos_from(key) ((key) & ref_id_bitmask)

static void dump_headers(FILE *fp, std::vector<std::string> &headers);
static void load_headers(FILE *fp, std::vector<std::string> &headers);
template <typename V>
static void dump_coordinates(FILE *fp, parlay::sequence<u8> &offsets, cqueue_t<V> &values);
template <typename V>
static void load_coordinates(FILE *fp, parlay::sequence<u8> &offsets, cqueue_t<V> &values);

template <typename K, typename V>
static u4 consolidate(cqueue_t<K> &q_keys, cqueue_t<V> &q_values, parlay::sequence<u8> &value_offsets, u4 block_sz);

////////////////////////////////////////////////////////////////////////////////

void j_index_t::init_query_buffers() {
    log_info("In j_index_t");
    hhs = new heavyhitter_ht_t<u4>[parlay::num_workers()];
    headers.emplace_back("*");
}

void j_index_t::add(std::string &name, parlay::slice<char *, char *> seq) {
    headers.push_back(name);
    auto frag_offset = frag_offsets.back();
    auto kmers = create_kmers(seq, k, sigma, encode_dna);
    const u4 stride = frag_len - frag_ovlp_len;
    u4 j = 0;
    for (u4 i = 0; i < kmers.size(); i += stride, j++) {
        size_t count = std::min(frag_len, (u4) kmers.size() - i);
        q_keys.push_back(kmers.data() + i, count);
        auto values = parlay::sequence<u4>(count, frag_offset + j);
        q_values.push_back(values.data(), count);
    }
    frag_offsets.push_back(frag_offset + j);
}

std::tuple<const char *, u4, float> j_index_t::search(parlay::slice<char *, char *> seq) {
    const auto i = parlay::worker_id();
    auto &hh = hhs[i];
    hh.reset();
    parlay::sequence<u4> keys = create_kmers_1t(seq, k, sigma, encode_dna);
    for (auto key : keys) {
        const auto &[vbegin, vend] = get(key);
        for (auto v = vbegin; v != vend; ++v) hh.insert(*v);
    }
    if (hh.top_key != -1) {
        float presence = (hh.top_count * 1.0) / keys.size();
        if (presence < presence_fraction) return {"*", 0, 0.0f};
        auto lb = lower_bound(frag_offsets.data(), 0, frag_offsets.size(), hh.top_key);
        auto &header = headers[lb - 1];
        u4 pos = (hh.top_key - lb) * (frag_len - frag_ovlp_len);
        return std::make_tuple(header.c_str(), pos, presence);
    } else return {"*", 0, 0.0f};
}

void j_index_t::build() {
    value_offsets.resize(n_keys+1);
    max_occ = consolidate(q_keys, q_values, value_offsets, sort_blocksz);
}

void j_index_t::dump(const std::string &filename) {
    FILE *fp = fopen(filename.c_str(), "a");
    expect(ftell(fp) == CONFIG_DUMP_SIZE);
    if (!fp) log_error("Cannot open file %s because %s.", filename.c_str(), strerror(errno));
    dump_headers(fp, headers);
    dump_coordinates(fp, value_offsets, q_values);
    fwrite(&max_occ, sizeof(max_occ), 1, fp);

    log_info("Dumping fragment offsets..");
    u4 buffer[3];
    buffer[0] = frag_offsets.size(), buffer[1] = frag_len, buffer[2] = frag_ovlp_len;
    fwrite(buffer, sizeof(u4), 3, fp);
    fwrite(frag_offsets.data(), sizeof(u4), buffer[0], fp);
    log_info("Done");
    fclose(fp);
}

void j_index_t::load(const std::string &filename) {
    FILE *fp = fopen(filename.c_str(), "r");
    if (!fp) log_error("Could not open %s because %s", filename.c_str(), strerror(errno));
    if (fseek(fp, CONFIG_DUMP_SIZE, SEEK_SET))
        log_error("Could not fseek because %s.", strerror(errno));
    load_headers(fp, headers);
    load_coordinates(fp, value_offsets, q_values);
    expect(fread(&max_occ, sizeof(max_occ), 1, fp) == 1);

    log_info("Loading fragment offsets..");
    u4 buffer[3];
    expect(fread(buffer, sizeof(u4), 3, fp) == 3);
    frag_offsets.resize(buffer[0]);
    frag_len = buffer[1], frag_ovlp_len = buffer[2];
    expect(fread(frag_offsets.data(), sizeof(u4), buffer[0], fp) == buffer[0]);
    log_info("Done");
    fclose(fp);
}

void c_index_t::init_query_buffers() {
    log_info("In c_index_t");
    hhs = new heavyhitter_ht_t<u8>[parlay::num_workers()];
}

void c_index_t::build() {
    value_offsets.resize(n_keys+1);
    max_occ = consolidate(q_keys, q_values, value_offsets, sort_blocksz);
}

void c_index_t::add(std::string &name, parlay::slice<char *, char *> seq) {
    const u4 id = headers.size();
    headers.push_back(name);

    auto kmers = create_kmers(seq, k, sigma, encode_dna);
    q_keys.push_back(kmers.data(), kmers.size());

    auto addresses = parlay::tabulate(kmers.size(), [&](size_t i) {
        return make_key_from(id, i);
    });

    q_values.push_back(addresses.data(), addresses.size());
}

std::tuple<const char *, u4, float> c_index_t::search(parlay::slice<char *, char *> seq) {
    const auto i = parlay::worker_id();
    auto &hh = hhs[i];
    hh.reset();
    parlay::sequence<u4> keys = create_kmers_1t(seq, k, sigma, encode_dna);
    for (u4 j = 0; j < keys.size(); ++j) {
        const auto &[vbegin, vend] = get(keys[j]);
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
    if (hh.top_key != -1) {
        float presence = (hh.top_count * 1.0) / keys.size();
        if (presence < presence_fraction) return {"*", 0, 0.0f};
        auto id = get_id_from(hh.top_key);
        auto &header = headers[get_id_from(hh.top_key)];
        return std::make_tuple(header.c_str(), get_pos_from(hh.top_key) * bandwidth, presence);
    } else return {"*", 0, 0.0f};
}

void c_index_t::dump(const std::string &filename) {
    FILE *fp = fopen(filename.c_str(), "a");
    expect(ftell(fp) == CONFIG_DUMP_SIZE);
    if (!fp) log_error("Cannot open file %s because %s.", filename.c_str(), strerror(errno));
    dump_headers(fp, headers);
    dump_coordinates(fp, value_offsets, q_values);
    fwrite(&max_occ, sizeof(max_occ), 1, fp);
    fclose(fp);
}

void c_index_t::load(const std::string &filename) {
    FILE *fp = fopen(filename.c_str(), "r");
    if (!fp) log_error("Could not open %s because %s", filename.c_str(), strerror(errno));
    if (fseek(fp, CONFIG_DUMP_SIZE, SEEK_SET))
        log_error("Could not fseek because %s.", strerror(errno));
    load_headers(fp, headers);
    load_coordinates(fp, value_offsets, q_values);
    expect(fread(&max_occ, sizeof(max_occ), 1, fp) == 1);
    fclose(fp);
}

template <typename K, typename V>
static void simple_sort_by_key(parlay::sequence<K> &keys, parlay::sequence<V> &values) {
    // Step 1: Create a vector of indices
    parlay::sequence<u4> indices(keys.size());
    for (int i = 0; i < keys.size(); ++i) {
        indices[i] = i; // Fill with the original indices
    }

    // Step 2: Sort the indices based on values in keys
    std::sort(indices.begin(), indices.end(), [&keys](int i1, int i2) {
        return keys[i1] < keys[i2];  // Compare based on the values of keys
    });

    // Step 3: Reorder keys and values according to the sorted indices
    parlay::sequence<K> keys_sorted(keys.size());
    parlay::sequence<V> values_sorted(keys.size());
    for (u4 i = 0; i < indices.size(); ++i) {
        keys_sorted[i] = keys[indices[i]];
        values_sorted[i] = values[indices[i]];
    }

    keys = std::move(keys_sorted);
    values = std::move(values_sorted);
}

template <typename T>
static parlay::sequence<u4> find_run_offsets(parlay::sequence<T> &values) {
    if (values.empty()) return {};
    parlay::sequence<u4> offsets(values.size()+1);
    offsets[0] = 0;
    u4 j = 1;
    for (size_t i = 1; i < values.size(); ++i) {
        if (values[i] != values[i - 1]) {
            offsets[j++] = i; // Mark the start of a new run
        }
    }
    offsets.resize(j+1);
    offsets.back() = values.size();
    return offsets;
}

static parlay::sequence<u4> get_minimizer_indices(parlay::sequence<u4> &keys, int w) {
    size_t n = keys.size();
    if (w > n) return {};
    auto midx = parlay::tabulate(n - w + 1, [&](size_t i){
        u4 min_key = -1, min_idx = -1;
        for (u4 j = i; j < i+w; ++j) {
            if (keys[j] < min_key) min_key = keys[j], min_idx = j;
        }
        return min_idx;
    });
    return parlay::unique(midx);
}

template <typename T>
static T calc_max_occ(parlay::sequence<T> &counts) {
    auto f_counts = parlay::filter(counts, [](T x) { return x != 0; });
    parlay::integer_sort_inplace(f_counts);
    auto occ = f_counts[f_counts.size() * 99/100];
    log_info("Min occ. = %lu", f_counts.front());
    log_info("Median occ. = %lu", f_counts[f_counts.size()/2]);
    log_info("Max occ. = %lu", f_counts.back());
    log_info("99%% = %lu", occ);
    return occ;
}

template <typename K, typename V>
static u4 consolidate(cqueue_t<K> &q_keys, cqueue_t<V> &q_values, parlay::sequence<u8> &value_offsets, u4 block_sz) {
    expect(q_keys.size() == q_values.size());
    log_info("Sorting %zd tuples..", q_keys.size());
    size_t bufsz = (sizeof(K) + sizeof(V)) * block_sz;
    auto buf = malloc(bufsz);
    cq_sort_by_key(q_keys, q_values, block_sz, buf);

    log_info("Counting unique keys..");
    std::vector<size_t> partition_sizes;
    cq_get_partitions(q_keys, block_sz, partition_sizes);
    verify(parlay::reduce(partition_sizes) == q_keys.size());
    auto d_keys_in = (K *) buf;
    for (auto np: partition_sizes) {
        q_keys.pop_front(d_keys_in, np);
        auto key_slice = parlay::slice(d_keys_in, d_keys_in + np);
        auto histogram = parlay::histogram_by_key(key_slice);
        parlay::sort_inplace(histogram);
        parlay::for_each(parlay::iota(histogram.size()), [&](size_t i) {
            auto key = histogram[i].first;
            value_offsets[key] = histogram[i].second;
        });
    }

    free(buf);
    MEMPOOL_SHRINK(K);
    MEMPOOL_SHRINK(V);
    PRINT_MEM_USAGE(K);
    PRINT_MEM_USAGE(V);
    auto occ99 = calc_max_occ(value_offsets);

    parlay::scan_inplace(value_offsets);
    PRINT_MEM_USAGE(K);
    PRINT_MEM_USAGE(V);
    return occ99;
}

static void dump_headers(FILE *fp, std::vector<std::string> &headers) {
    size_t n = headers.size();
    log_info("Dumping %zd refnames..", n);
    fwrite(&n, sizeof(n), 1, fp);
    std::vector<u2> refnamelens(n);
    for (int i = 0; i < n; ++i) refnamelens[i] = headers[i].size();
    fwrite(refnamelens.data(), sizeof(u2), n, fp);
    for (const auto& refname : headers) {
        n = refname.size();
        fwrite(refname.c_str(), 1, n, fp);
    }
    log_info("Done.");
}

static void load_headers(FILE *fp, std::vector<std::string> &headers) {
    size_t n;
    expect(fread(&n, sizeof(n), 1, fp) == 1);
    log_info("Loading %zd refnames..", n);
    std::vector<u2> refnamelens(n);
    expect(fread(refnamelens.data(), sizeof(u2), n, fp) == n);

    char buffer[4096];
    for (int i = 0; i < n; ++i) {
        size_t l = refnamelens[i];
        expect(fread(buffer, 1, l, fp) == 1);
        headers.emplace_back(buffer, l);
    }
    log_info("Done.");
}

template <typename V>
static void dump_coordinates(FILE *fp, parlay::sequence<u8> &offsets, cqueue_t<V> &values) {
    size_t n_keys = offsets.size(), n_values = values.size();
    log_info("Dumping counts..");
    fwrite(&n_keys, sizeof(n_keys), 1, fp);
    fwrite(offsets.data(), sizeof(u8), offsets.size(), fp);

    log_info("Dumping values..");
    values.dump(fp);

    log_info("Done.");
}

template <typename V>
static void load_coordinates(FILE *fp, parlay::sequence<u8> &offsets, cqueue_t<V> &values) {
    size_t n_keys, n_values;
    log_info("Loading counts..");
    expect(fread(&n_keys, sizeof(n_keys), 1, fp) == 1);
    offsets.resize(n_keys);
    expect(fread(offsets.data(), sizeof(u8), offsets.size(), fp) == offsets.size());

    log_info("Loading values..");
    values.load(fp);

    log_info("Done.");
}
