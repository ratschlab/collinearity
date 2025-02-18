//
// Created by Sayan Goswami on 18.02.2025.
//
#include "collinearity.h"

jindex_t::jindex_t(int k, int sigma): k(k), sigma(sigma) {}

void jindex_t::init_query_buffers() {
    hhs = new heavyhitter_ht_t<u4>[parlay::num_workers()];
    headers.emplace_back("*");
}

std::tuple<const char*, u4, float> jindex_t::get(std::string &seq) {
    const size_t n = seq.size(), n_keys = n - k + 1;
    parlay::sequence<u4> keys = create_kmers_1t(seq, k, sigma, encode_dna);
    expect(keys.size() == n_keys);

    const auto i = parlay::worker_id();
    auto &hh = hhs[i];
    hh.reset();

    for (u4 j = 0; j < n_keys; ++j) {
        const auto &[vbegin, vend] = get(keys[j]);
        for (auto v = vbegin; v != vend; ++v) {
            hh.insert(*v);
        }
    }

    // in this case, presence_fraction is the frac. of kmers that support an intercept
    float presence = (hh.top_count * 1.0) / n_keys;
    if (presence >= presence_fraction) {
        auto lb = lower_bound(frag_offsets.data(), 0, frag_offsets.size(), hh.top_key);
        auto &header = headers[lb];
        lb = lb? lb - frag_offsets[lb-1] : lb;
        u4 pos = (hh.top_key - lb) * (frag_len - frag_ovlp_len);
        return std::make_tuple(header.c_str(), pos, presence);
    }
    else return std::make_tuple(headers.back().c_str(), -1, 0.0f);
}

void jindex_t::dump(const std::string &basename) {
    throw "Not implemented";
}

void jindex_t::load(const std::string &basename) {
    throw "Not implemented";
}

std::pair<u4 *, u4 *> jindex_t::get(u4 key) {
    u4 shard_id = SHARD(key);
    auto &shard = shards[shard_id];
    if (shard.tuples.contains(SKEY(key))) {
        auto val = shard.tuples[SKEY(key)];
        auto offset = val >> 32, count = val & 0xffffffff;
//            if (count >= max_allowed_occ) return {nullptr, nullptr};
        return {shard.values.data() + offset, shard.values.data() + offset + count};
    } else return {nullptr, nullptr};
}
