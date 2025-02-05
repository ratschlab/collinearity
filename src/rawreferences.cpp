//
// Created by Sayan Goswami on 23.01.2025.
//

#include "collinearity.h"

using namespace std;
using namespace parlay;

pair<int, vector<double>> load_pore_model(string &poremodel_file) {
    vector<string> kmers;
    vector<double> levels;
    ifstream file(poremodel_file);
    if (!file.is_open()) log_error("Could not open file %s because %s.", poremodel_file.c_str(), strerror(errno));

    string line;
    // Read the header line
    if (!getline(file, line)) log_error("File is empty or missing header.");

    // Parse each row
    while (getline(file, line)) {
        istringstream iss(line);
        string token, kmer;
        double level;

        // Parse kmer
        if (!getline(iss, kmer, '\t')) log_error("Error parsing kmer.");

        // Parse level_mean
        if (!getline(iss, token, '\t') || !(istringstream(token) >> level))
            log_error("Error parsing level_mean.");

        // Add the row to the vector
        kmers.push_back(kmer);
        levels.push_back(level);
    }

    int k = kmers[0].length();
    log_info("Poremodel K = %d.", k);

    // Encode kmers
    vector<u4> kmerints;
    for (const auto& kmer : kmers) kmerints.push_back(encode_kmer(kmer, k, 4, encode_dna));

    // Check if all k-mers are represented
    int exptd_num_kmers = 1 << (k << 1);
    int num_kmers = kmerints.size();
    expect(num_kmers == exptd_num_kmers);
    expect(*min_element(kmerints.begin(), kmerints.end()) == 0);
    expect(*max_element(kmerints.begin(), kmerints.end()) == exptd_num_kmers - 1);

    vector<double> levels1(levels.size());
    for (int i = 0; i < num_kmers; ++i) levels1[kmerints[i]] = levels[i];

    sort(kmerints.begin(), kmerints.end());
    vector<int> diffs(num_kmers - 1);
    adjacent_difference(kmerints.begin(), kmerints.end(), diffs.begin());
    if (all_of(diffs.begin() + 1, diffs.end(), [](int diff) { return diff == 1; })) log_info("Done");
    else log_error("kmers missing from pore model");

    // normalize levels
    double sum = accumulate(levels1.begin(), levels1.end(), 0.0);
    double mean = sum / num_kmers;
    vector<double> diff(num_kmers);
    transform(levels1.begin(), levels1.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = sqrt(sq_sum / num_kmers);
    log_info("Mean = %f, StdDev = %f.", mean, stdev);

    transform(levels1.begin(), levels1.end(), levels.begin(),
              [mean, stdev](double val) { return (val - mean) / stdev; });

    return {k, levels};
}

template <typename T>
parlay::sequence<double> sequence2squiggles(const T &sequence, const int k, const vector<double> &levels) {
    auto kmers = create_kmers(sequence, k, 4, encode_dna);
    return parlay::map(iota(kmers.size()), [&](size_t i) {
        return levels[kmers[i]];
    });
}

template parlay::sequence<double> sequence2squiggles<std::string>(const std::string &sequence, const int k, const vector<double> &levels);
template parlay::sequence<double> sequence2squiggles<parlay::sequence<char>>(const parlay::sequence<char> &sequence, const int k, const vector<double> &levels);
