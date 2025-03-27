//
// Created by Sayan Goswami on 18.02.2025.
//
#include "collinearity.h"
#include "index.h"


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
template <typename K, typename V>
static void dump_hashmap(FILE *fp, shard_t<K,V> *shards);
template <typename K, typename V>
static void load_hashmap(FILE *fp, shard_t<K,V> *shards);
template <typename K, typename V>
static void dump_values(FILE *fp, shard_t<K,V> *shards);
template <typename K, typename V>
static void load_values(FILE *fp, shard_t<K,V> *shards);

////////////////////////////////////////////////////////////////////////////////

template<typename K, typename V>
std::pair<V *, V *> get(shard_t<K,V> *shards, K key) {
    u4 shard_id = SHARD(key);
    auto &shard = shards[shard_id];
    if (shard.tuples.contains(SKEY(key))) {
        auto val = shard.tuples[SKEY(key)];
        auto offset = val >> 32, count = val & 0xffffffff;
//            if (count >= max_allowed_occ) return {nullptr, nullptr};
        return {shard.values.data() + offset, shard.values.data() + offset + count};
    } else return {nullptr, nullptr};
}

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
        const auto &[vbegin, vend] = get(shards, key);
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
    max_occ = consolidate(q_keys, q_values, shards);
}

void j_index_t::dump(const std::string &filename) {
    FILE *fp = fopen(filename.c_str(), "a");
    expect(ftell(fp) == CONFIG_DUMP_SIZE);
    if (!fp) log_error("Cannot open file %s because %s.", filename.c_str(), strerror(errno));
    dump_headers(fp, headers);
    dump_hashmap(fp, shards);
    dump_values(fp, shards);
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
    load_hashmap(fp, shards);
    load_values(fp, shards);
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
    max_occ = consolidate(q_keys, q_values, shards);
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
        const auto &[vbegin, vend] = get(shards, keys[j]);
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
    dump_hashmap(fp, shards);
    dump_values(fp, shards);
    fwrite(&max_occ, sizeof(max_occ), 1, fp);

    fclose(fp);
}

void c_index_t::load(const std::string &filename) {
    FILE *fp = fopen(filename.c_str(), "r");
    if (!fp) log_error("Could not open %s because %s", filename.c_str(), strerror(errno));
    if (fseek(fp, CONFIG_DUMP_SIZE, SEEK_SET))
        log_error("Could not fseek because %s.", strerror(errno));
    load_headers(fp, headers);
    load_hashmap(fp, shards);
    load_values(fp, shards);
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



template <typename K, typename V>
static void put_in_shard(shard_t<K, V> &shard, parlay::sequence<K> &all_keys, parlay::sequence<V> &all_values, parlay::sequence<u8> &my_indices) {
    u4 p = my_indices.size();
    parlay::sequence<K> my_keys(p);
    parlay::sequence<V> my_values(p);
    for (u4 i = 0; i < p; ++i) {
        my_keys[i] = SKEY(all_keys[my_indices[i]]);
        my_values[i] = all_values[my_indices[i]];
    }
    simple_sort_by_key(my_keys, my_values);
    auto offsets = find_run_offsets(my_keys);

    size_t old_n_vals = shard.values.size();
    shard.values.append(my_values);

    for (int i = 0; i < offsets.size()-1; ++i) {
        auto key = my_keys[offsets[i]];
        auto n_vals = offsets[i+1] - offsets[i];
        auto offset = old_n_vals + offsets[i];
        shard.tuples[key] = ((offset << 32) | n_vals);
    }
}

template <typename K, typename V>
static u4 calc_max_occ(shard_t<K,V> *shards) {
    u4 max_allowed_occ = -1;
    parlay::sequence<u4> n_unique_keys_per_shard(N_SHARDS);
    parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
        n_unique_keys_per_shard[i] = shards[i].tuples.size();
    });
    u4 total_uniq_keys = parlay::reduce(n_unique_keys_per_shard);
    log_info("# unique kmers = %u", total_uniq_keys);
    if (total_uniq_keys) {
        parlay::sequence<u4> occ(total_uniq_keys);
        auto [offset, total] = parlay::scan(n_unique_keys_per_shard);
        expect(total == total_uniq_keys);

        parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i) {
            auto my_off = offset[i];
            int j = 0;
            for (auto &[key, value]: shards[i].tuples) {
                occ[my_off + j] = value & 0xffffffff;
                j++;
            }
        });
        parlay::sort_inplace(occ);
        log_info("Min occ. = %u", occ[0]);
        log_info("Median occ. = %u", occ[total_uniq_keys / 2]);
        log_info("Max occ. = %u", occ[total_uniq_keys - 1]);
        max_allowed_occ = occ[total_uniq_keys * .99];
        log_info("99%% = %u", max_allowed_occ);
    } else max_allowed_occ = 1<<31;
    return max_allowed_occ;
}

template <typename K, typename V>
static u4 consolidate(cqueue_t<K> &q_keys, cqueue_t<V> &q_values, shard_t<K, V> *shards) {
    expect(q_keys.size() == q_values.size());
    log_info("Sorting %zd tuples..", q_keys.size());
    size_t bufsz = (sizeof(K) + sizeof(V)) * BLOCK_SZ;
    auto buf = malloc(bufsz);
    cq_sort_by_key(q_keys, q_values, BLOCK_SZ, buf);

    std::vector<size_t> partition_sizes;
    cq_get_partitions(q_keys, BLOCK_SZ, partition_sizes);
    verify(parlay::reduce(partition_sizes) == q_keys.size());

    parlay::sequence<K> tmp_keys(BLOCK_SZ);
    parlay::sequence<V> tmp_values(BLOCK_SZ);

    for (auto np : partition_sizes) {
        expect(np);
        tmp_keys.resize(np), tmp_values.resize(np);
        q_keys.pop_front(tmp_keys.data(), np);
        q_values.pop_front(tmp_values.data(), np);

        auto tuples = parlay::delayed_tabulate(np, [&](size_t i){
            return std::make_pair(tmp_keys[i], i);
        });
        auto hash = [](const u4 &x) { return SHARD(x); };
        auto equal = [](const u4& a, const u4& b) { return SHARD(a) == SHARD(b); };
        auto grouped = parlay::group_by_key(tuples, hash, equal);

        parlay::for_each(parlay::iota(grouped.size()), [&](size_t i){
            auto shard_id = (SHARD(grouped[i].first));
            put_in_shard(shards[shard_id], tmp_keys, tmp_values, grouped[i].second);
        });
    }
    return calc_max_occ(shards);
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

template <typename K, typename V>
static void dump_hashmap(FILE *fp, shard_t<K,V> *shards) {
    log_info("Collecting key-value pairs from shards..");
    parlay::sequence<u4> n_unique_keys_per_shard(N_SHARDS);
    parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
        n_unique_keys_per_shard[i] = shards[i].tuples.size();
    });
    u4 total_uniq_keys = parlay::reduce(n_unique_keys_per_shard);
    parlay::sequence<K> keys(total_uniq_keys);
    parlay::sequence<u8> vpartitions(total_uniq_keys);
    auto [offset, total] = parlay::scan(n_unique_keys_per_shard);
    expect(total == total_uniq_keys);
    offset.push_back(total_uniq_keys);

    parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
        auto my_off = offset[i];
        int j = 0;
        for (auto &[key, vpartition]: shards[i].tuples) {
            keys[my_off + j] = key;
            vpartitions[my_off + j] = vpartition;
            j++;
        }
        expect(j == offset[i+1] - offset[i]);
    });

    log_info("Dumping keys..");
    expect(offset.size() == N_SHARDS + 1);
    fwrite(&total_uniq_keys, sizeof(total_uniq_keys), 1, fp);
    fwrite(offset.data(), sizeof(offset[0]), offset.size(), fp);
    fwrite(keys.data(), sizeof(keys[0]), keys.size(), fp);
    fwrite(vpartitions.data(), sizeof(vpartitions[0]), vpartitions.size(), fp);
    log_info("Done");
}

template <typename K, typename V>
static void load_hashmap(FILE *fp, shard_t<K,V> *shards) {
    log_info("Loading keys..");
    u4 total_uniq_keys;
    parlay::sequence<u4> offset(N_SHARDS+1);
    expect(fread(&total_uniq_keys, sizeof(total_uniq_keys), 1, fp) == 1);
    expect(fread(offset.data(), sizeof(offset[0]), offset.size(), fp) == offset.size());

    parlay::sequence<K> keys(total_uniq_keys);
    parlay::sequence<u8> vpartitions(total_uniq_keys);
    expect(fread(keys.data(), sizeof(keys[0]), keys.size(), fp) == keys.size());
    expect(fread(vpartitions.data(), sizeof(vpartitions[0]), vpartitions.size(), fp) == vpartitions.size());

    log_info("Adding keys into hash table");
    parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
        auto my_off = offset[i], count = offset[i+1] - offset[i];
        for (int j = 0; j < count; ++j)
            shards[i].tuples[keys[my_off + j]] = vpartitions[my_off + j];
    });
    log_info("Done");
}

template <typename K, typename V>
static void dump_values(FILE *fp, shard_t<K,V> *shards) {
    log_info("Dumping values..");
    for (u4 i = 0; i < N_SHARDS; ++i) {
        u4 n_values = shards[i].values.size();
        fwrite(&n_values, sizeof(n_values), 1, fp);
        fwrite(shards[i].values.data(), sizeof(shards[i].values[0]), shards[i].values.size(), fp);
    }
    log_info("Done.");
}

template <typename K, typename V>
static void load_values(FILE *fp, shard_t<K,V> *shards) {
    log_info("Loading values..");
    for (u4 i = 0; i < N_SHARDS; ++i) {
        u4 n_values = shards[i].values.size();
        expect(fread(&n_values, sizeof(n_values), 1, fp) == 1);
        shards[i].values.resize(n_values);
        expect(fread(shards[i].values.data(), sizeof(shards[i].values[0]), shards[i].values.size(), fp) == shards[i].values.size());
    }
    log_info("done.");
}