//
// Created by Sayan Goswami on 23.01.2025.
//

#include "collinearity.h"

using namespace std;
//using namespace parlay;

static pair<parlay::sequence<double>, parlay::sequence<double>> compute_prefix_prefixsq(const parlay::sequence<double> &signal) {
    const auto s_len = signal.size();
    parlay::sequence<double> prefix_sum(s_len+1, 0), prefix_sum_square(s_len+1, 0);
    inclusive_scan(signal.begin(), signal.end(), prefix_sum.begin()+1);
    transform_inclusive_scan(signal.begin(), signal.end(), prefix_sum_square.begin()+1, plus<>(),
                             [](double x) { return x * x; });
    return {prefix_sum, prefix_sum_square};
}

static parlay::sequence<double> compute_tstat(parlay::sequence<double> &prefix_sum, parlay::sequence<double> &prefix_sum_square, int w) {
    // Compute windowed tstat
    const auto s_len = prefix_sum.size() - 1;  // length of original signal
    double eta = numeric_limits<double>::epsilon();  // Equivalent of FLT_MIN
    parlay::sequence<double> tstat(s_len + 1, 0.0);

    // Check for conditions
    if (s_len < 2 * w || w < 2) {
        return tstat;
    }

    for (int i = w; i <= s_len - w; ++i) {
        double sum1 = prefix_sum[i];
        double sumsq1 = prefix_sum_square[i];
        if (i > w) {
            sum1 -= prefix_sum[i - w];
            sumsq1 -= prefix_sum_square[i - w];
        }

        double sum2 = prefix_sum[i + w] - prefix_sum[i];
        double sumsq2 = prefix_sum_square[i + w] - prefix_sum_square[i];

        double mean1 = sum1 / w;
        double mean2 = sum2 / w;

        double combined_var = (sumsq1 / w - mean1 * mean1 + sumsq2 / w - mean2 * mean2) / w;

        // Ensure combined_var is not too small
        combined_var = max(combined_var, eta);

        // Compute t-statistic
        double delta_mean = mean2 - mean1;
        tstat[i] = abs(delta_mean) / sqrt(combined_var);
    }

    // Set last `w_len` elements to 0 (fudging boundaries)
    fill(tstat.begin() + s_len - w + 1, tstat.end(), 0.0f);

    return tstat;
}

parlay::sequence<size_t> tstat_segmenter_t::segment_signal(const parlay::sequence<double> &signal) {
    const float peak_height = 0.4f;
    parlay::sequence<size_t> peaks;
    const int n_detectors = detectors.size();
    const auto s_len = signal.size();
    auto [prefix_sum, prefix_sum_sq] = compute_prefix_prefixsq(signal);
    for(auto &detector: detectors) {
        detector.reset();
        detector.tstat = compute_tstat(prefix_sum, prefix_sum_sq, detector.WINDOW_LENGTH);
    }

    for (uint32_t i = 0; i < s_len; i++) {
        for (uint32_t k = 0; k < n_detectors; k++) {
            auto &detector = detectors[k];
            if (detector.masked_to >= i) continue;

            double current_value = detector.tstat[i];
            // double adaptive_peak_height = calculate_adaptive_peak_height(prefix_sum, prefix_sum_square, i, detector->window_length, peak_height);

            if (detector.peak_pos == detector.DEF_PEAK_POS) {
                // CASE 1: We've not yet recorded any maximum
                if (current_value < detector.peak_value) { // A deeper minimum:
                    detector.peak_value = current_value;
                } else if (current_value - detector.peak_value > peak_height) {
                    // ...or a qualifying maximum:
                    detector.peak_value = current_value;
                    detector.peak_pos = i;
                    // otherwise, wait to rise high enough to be considered a peak
                }
            } else {
                // CASE 2: In an existing peak, waiting to see if it is good
                if (current_value > detector.peak_value) {
                    // Update the peak
                    detector.peak_value = current_value;
                    detector.peak_pos = i;
                }
                // Tell other detectors no need to check for a peak until a certain point
                if (detector.peak_value > detector.THRESHOLD) {
                    for(int n_d = k+1; n_d < n_detectors; n_d++){
                        detectors[n_d].masked_to = detector.peak_pos + detectors[0].WINDOW_LENGTH;
                        detectors[n_d].peak_pos = detectors[n_d].DEF_PEAK_POS;
                        detectors[n_d].peak_value = detectors[n_d].DEF_PEAK_VAL;
                        detectors[n_d].valid_peak = 0;
                    }
                }
                // There is a good peak
                if (detector.peak_value - current_value > peak_height &&
                    detector.peak_value > detector.THRESHOLD) {
                    detector.valid_peak = 1;
                }
                // Check if we are now further away from the current peak
                if (detector.valid_peak && (i - detector.peak_pos) > detector.WINDOW_LENGTH / 2) {
                    peaks.push_back(detector.peak_pos);
                    detector.peak_pos = detector.DEF_PEAK_POS;
                    detector.peak_value = current_value;
                    detector.valid_peak = 0;
                }
            }
        }
    }

    return peaks;
}