//
// Created by Sayan Goswami on 09.12.2024.
//

#ifndef COLLINEARITY_INDEX_H
#define COLLINEARITY_INDEX_H

#include "prelude.h"
#include "parlay_utils.h"
#include "cqueue.h"
#include "cqutils.h"
#include "hash_table8.hpp"
#include "utils.h"

#ifdef NDEBUG
#define SANITY_CHECKS 0
#else
#define SANITY_CHECKS 1
#endif

// I am assuming that the number of sequences won't exceed 2^20 (1M)
// and the longest sequence won't be longer than 2^40 (1T)
// if that's not the case, change the following line accordingly
#define ref_id_nbits 20
#define ref_len_nbits (64 - ref_id_nbits)
#define ref_id_bitmask (((1ULL) << ref_len_nbits) - 1)
#define make_key_from(id, pos) (((u8)(id)) << (ref_len_nbits) | (pos))
#define get_id_from(key) ((key) >> ref_len_nbits)
#define get_pos_from(key) ((key) & ref_id_bitmask)

#define N_SHARD_BITS 7
#define N_SHARDS (1<<N_SHARD_BITS)
#define SHARD(x) ((x) & (N_SHARDS-1))
#define SKEY(x) ((x) >> N_SHARD_BITS)

// the minimum fraction of kmers required for majority vote
const float presence_fraction = .1f;

// this is the bandwidth, I will later set this from the config
const int bandwidth = 15;

template <typename T>
struct [[maybe_unused]] heavyhitter_ht_t {
    emhash8::HashMap<T,u4> counts;
    T top_key = -1;
    u4 top_count = 0;
    void insert(const T key) {
        u4 count = counts[key]++;
        if (count > top_count)
            top_count = count, top_key = key;
    }
    void reset() { counts.clear(), top_key=-1, top_count = 0; }
};

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
struct shard_t {
    emhash8::HashMap<K,u8> tuples;
    parlay::sequence<V> values;
};

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
    log_info("Sorting..");
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

class jindex_t {
    shard_t<u4,u4> shards[N_SHARDS];
    cqueue_t<u4> q_keys;
    cqueue_t<u4> q_values;
    heavyhitter_ht_t<u4> *hhs = nullptr;
    const u4 k, sigma;
    const u4 frag_len = 180, frag_ovlp_len = 120;
    std::vector<std::string> headers;
    std::vector<u4> frag_offsets = {0};
    u4 max_occ = -1;

    std::pair<u4*, u4*> get(u4 key);

public:
    explicit jindex_t(int k, int sigma);
    void init_query_buffers();
    void build() { max_occ = consolidate(q_keys, q_values, shards); }
    std::tuple<const char*, u4, float> get(std::string &seq);
    void dump(const std::string& basename);
    void load(const std::string& basename);

    template <typename Sequence>
    void put(std::string &name, Sequence &seq) {
        headers.push_back(name);
        auto frag_offset = frag_offsets.back();

        auto kmers = create_kmers(seq, k, sigma, encode_dna);
        const u4 stride = frag_len - frag_ovlp_len;
        u4 j = 0;
        for (u4 i = 0; i < kmers.size(); i += stride, j++) {
            size_t count =  std::min(frag_len, (u4)kmers.size() - i);
            q_keys.push_back(kmers.data() + i, count);
            auto values = parlay::sequence<u4>(count, frag_offset + j);
            q_values.push_back(values.data(), count);
        }
        frag_offsets.push_back(frag_offset + j);
    }
};

struct sindex_t {
    shard_t<u4,u8> shards[N_SHARDS];
    cqueue_t<u4> q_keys;
    cqueue_t<u8> q_values;
    std::vector<std::string> headers;
    u4 max_allowed_occ;
    heavyhitter_ht_t<u8> *hhs = nullptr;
    const bool minimize = true;

    void add(std::string &name, parlay::sequence<u4> &keys, parlay::sequence<u8> &addresses) {
        headers.push_back(name);
        q_keys.push_back(keys.data(), keys.size());
        q_values.push_back(addresses.data(), addresses.size());
    }

    void build() {
        max_allowed_occ = consolidate(q_keys, q_values, shards);
    }

    void init_query_buffers() {hhs = new heavyhitter_ht_t<u8>[parlay::num_workers()];}

    std::tuple<u4, u4, float> search(const u4 *keys, const size_t n_keys) {
        const auto i = parlay::worker_id();
        auto &hh = hhs[i];
        hh.reset();

        for (u4 j = 0; j < n_keys; ++j) {
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

        // in this case, presence_fraction is the frac. of kmers that support an intercept
        float presence = (hh.top_count * 1.0) / n_keys;
        if (presence >= presence_fraction)
            return std::make_tuple(get_id_from(hh.top_key), get_pos_from(hh.top_key) * bandwidth, presence);
        else return std::make_tuple(-1, -1, 0.0f);
    }

    std::pair<u8*, u8*> get(u4 key) {
        u4 shard_id = SHARD(key);
        auto &shard = shards[shard_id];
        if (shard.tuples.contains(SKEY(key))) {
            auto val = shard.tuples[SKEY(key)];
            auto offset = val >> 32, count = val & 0xffffffff;
//            if (count >= max_allowed_occ) return {nullptr, nullptr};
            return {shard.values.data() + offset, shard.values.data() + offset + count};
        } else return {nullptr, nullptr};
    }

    void dump(const std::string& basename) {
        log_info("Dumping index.");
        std::string filename = basename + ".cidx";
        FILE *fp = fopen(filename.c_str(), "wb");
        if (!fp) log_error("Could not open %s because %s", filename.c_str(), strerror(errno));

        size_t n = headers.size();
        log_info("Dumping %zd refnames..", n);
        fwrite(&n, sizeof(n), 1, fp);
        std::vector<u2> refnamelens(n);
        for (int i = 0; i < n; ++i) refnamelens[i] = headers[i].size();
        fwrite(refnamelens.data(), sizeof(u2), n, fp);
        for (const auto& refname : headers) {
            n = refname.size();
            fwrite(refname.c_str(), n, 1, fp);
        }

        log_info("Collecting key-value pairs from shards..");
        parlay::sequence<u4> n_unique_keys_per_shard(N_SHARDS);
        parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
            n_unique_keys_per_shard[i] = shards[i].tuples.size();
        });
        u4 total_uniq_keys = parlay::reduce(n_unique_keys_per_shard);
        parlay::sequence<u4> keys(total_uniq_keys);
        parlay::sequence<u8> values(total_uniq_keys);
        auto [offset, total] = parlay::scan(n_unique_keys_per_shard);
        expect(total == total_uniq_keys);
        offset.push_back(total_uniq_keys);

        parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
            auto my_off = offset[i];
            int j = 0;
            for (auto &[key, value]: shards[i].tuples) {
                keys[my_off + j] = key;
                values[my_off + j] = value;
                j++;
            }
            expect(j == offset[i+1] - offset[i]);
        });

        log_info("Dumping keys..");
        expect(offset.size() == N_SHARDS + 1);
        fwrite(&total_uniq_keys, sizeof(total_uniq_keys), 1, fp);
        fwrite(offset.data(), sizeof(offset[0]), offset.size(), fp);
        fwrite(keys.data(), sizeof(keys[0]), keys.size(), fp);
        fwrite(values.data(), sizeof(values[0]), values.size(), fp);

        log_info("Dumping values..");
        for (auto & shard : shards) {
            u4 n_values = shard.values.size();
            fwrite(&n_values, sizeof(n_values), 1, fp);
            fwrite(shard.values.data(), sizeof(shard.values[0]), shard.values.size(), fp);
        }

        fwrite(&max_allowed_occ, sizeof(max_allowed_occ), 1, fp);

        fclose(fp);
        log_info("Done.");
    }

    void load(const std::string& basename) {
        log_info("Loading index.");
        std::string filename = basename + ".cidx";
        FILE *fp = fopen(filename.c_str(), "rb");
        if (!fp) log_error("Could not open %s because %s", filename.c_str(), strerror(errno));

        size_t n;
        expect(fread(&n, sizeof(n), 1, fp) == 1);
        log_info("Loading %zd refnames..", n);
        std::vector<u2> refnamelens(n);
        expect(fread(refnamelens.data(), sizeof(u2), n, fp) == n);

        char buffer[4096];
        for (int i = 0; i < n; ++i) {
            size_t l = refnamelens[i];
            expect(fread(buffer, l, 1, fp) == 1);
            headers.emplace_back(buffer, l);
        }

        log_info("loading keys..");
        u4 total_uniq_keys;
        parlay::sequence<u4> offset(N_SHARDS+1);
        expect(fread(&total_uniq_keys, sizeof(total_uniq_keys), 1, fp) == 1);
        expect(fread(offset.data(), sizeof(offset[0]), offset.size(), fp) == offset.size());

        parlay::sequence<u4> keys(total_uniq_keys);
        parlay::sequence<u8> values(total_uniq_keys);
        expect(fread(keys.data(), sizeof(keys[0]), keys.size(), fp) == keys.size());
        expect(fread(values.data(), sizeof(values[0]), values.size(), fp) == values.size());

        log_info("adding keys into hash table");
        parlay::for_each(parlay::iota(N_SHARDS), [&](size_t i){
            auto my_off = offset[i], count = offset[i+1] - offset[i];
            for (int j = 0; j < count; ++j)
                shards[i].tuples[keys[my_off + j]] = values[my_off + j];
        });

        log_info("loading values..");
        for (auto & shard : shards) {
            u4 n_values = shard.values.size();
            expect(fread(&n_values, sizeof(n_values), 1, fp) == 1);
            shard.values.resize(n_values);
            expect(fread(shard.values.data(), sizeof(shard.values[0]), shard.values.size(), fp) == shard.values.size());
        }

        expect(fread(&max_allowed_occ, sizeof(max_allowed_occ), 1, fp) == 1);

        fclose(fp);
        log_info("Done.");
    }
};

static inline parlay::sequence<u8> create_addresses(u8 id, u8 num_kmers) {
    return parlay::tabulate(num_kmers, [&](size_t i) {
        return make_key_from(id, i);
    });
}

#endif //COLLINEARITY_INDEX_H
