#include "collinearity.h"

static void print_idx_info(index_t *idx);
std::pair<index_t*, std::vector<std::string>> load_index(char *filename);
void dump_index(index_t *idx, std::vector<std::string> &refnames, const char *basename);

int main(int argc, char *argv[]) {
    info("Hello!");
    const char *ref = "/scratch/Zymo/Refs-fwd-rev.fasta";
    const char *qry = "/scratch/Zymo/basecalled/fast/Sigs.fasta"; //"/scratch/Zymo/reads-tiny.fasta";

    if (str_endswith(argv[1], ".cidx")) {
        auto [idx, refnames] = load_index(argv[1]);
        if (SANITY_CHECKS) print_idx_info(idx);
        query(argv[2], KMER_LENGTH, SIGMA, BATCH_SZ, idx, refnames);
    } else {
        auto [idx, refnames] = process_fasta(argv[1], KMER_LENGTH, SIGMA);
        dump_index(idx, refnames, argv[1]);
        if (SANITY_CHECKS) print_idx_info(idx);
        if (argc > 2) query(argv[2], KMER_LENGTH, SIGMA, BATCH_SZ, idx, refnames);
    }
    info("Bye!");
    return 0;
}

void dump_index(index_t *idx, std::vector<std::string> &refnames, const char *basename) {
    info("Dumping index.");
    std::string filename = std::string(basename) + std::string(".cidx");
    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp) error("Could not open %s because %s", filename.c_str(), strerror(errno));

    size_t n = refnames.size();
    info("Dumping %zd refnames..", n);
    fwrite(&n, sizeof(n), 1, fp);
    std::vector<u2> refnamelens(n);
    for (int i = 0; i < n; ++i) refnamelens[i] = refnames[i].size();
    fwrite(refnamelens.data(), sizeof(u2), n, fp);
    for (const auto& refname : refnames) {
        n = refname.size();
        fwrite(refname.c_str(), n, 1, fp);
    }

    n = idx->counts.size();
    info("Dumping %zd counts and offsets..", n);
    fwrite(&n, sizeof(n), 1, fp);
    fwrite(idx->counts.data(), sizeof(u4), n, fp);
    fwrite(idx->offsets.data(), sizeof(u8), n, fp);

    n = idx->values.size();
    info("Dumping %zd values..", n);
    fwrite(&n, sizeof(n), 1, fp);
    fwrite(idx->values.data(), sizeof(u8), n, fp);
    fclose(fp);
    info("Done.");
}

std::pair<index_t*, std::vector<std::string>> load_index(char *filename) {
    info("Loading index.");
    FILE *fp = fopen(filename, "rb");
    if (!fp) error("Could not open %s because %s", filename, strerror(errno));

    size_t n;
    fread(&n, sizeof(n), 1, fp);
    info("Loading %zd refnames..", n);
    std::vector<u2> refnamelens(n);
    fread(refnamelens.data(), sizeof(u2), n, fp);

    std::vector<std::string> refnames;
    char buffer[4096];
    for (int i = 0; i < n; ++i) {
        size_t l = refnamelens[i];
        fread(buffer, l, 1, fp);
        refnames.emplace_back(buffer, l);
    }

    auto *idx = new index_t;
    fread(&n, sizeof(n), 1, fp);
    info("Loading %zd counts and offsets..", n);
    idx->counts.resize(n);
    idx->offsets.resize(n);
    fread(idx->counts.data(), sizeof(u4), n, fp);
    fread(idx->offsets.data(), sizeof(u8), n, fp);

    fread(&n, sizeof(n), 1, fp);
    info("Loading %zd values..", n);
    idx->values.resize(n);
    fread(idx->values.data(), sizeof(u8), n, fp);

    fclose(fp);
    return std::make_pair(idx, std::move(refnames));
}

static void print_idx_info(index_t *idx) {
    auto counts_copy = parlay::integer_sort(idx->counts);
    auto start = lower_bound<u4>(counts_copy.data(), 0, counts_copy.size(), 1);
    size_t nk = counts_copy.size() - start;
    info("# kmers = %zd", nk);
    info("Min occ. = %u", counts_copy[start]? counts_copy[start] : counts_copy[start+1]);
    info("Median occ. = %u", counts_copy[start + nk / 2]);
    info("Max occ. = %u", counts_copy.back());
}