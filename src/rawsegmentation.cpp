//
// Created by Sayan Goswami on 23.01.2025.
//

#include "collinearity.h"

static double calculate_mean_of_filtered_segment(double *begin, double *end) {
    // Calculate median and IQR
    std::sort(begin, end);  // Sort the segment
    size_t segment_length = end - begin;
    double q1 = begin[segment_length / 4];
    double q3 = begin[3 * segment_length / 4];
    double iqr = q3 - q1;
    double lower_bound = q1 - iqr;
    double upper_bound = q3 + iqr;

    double sum_values = 0.0;
    int count = 0;
    for (auto it = begin; it != end; ++it) {
        auto value = *it;
        if (lower_bound <= value && value <= upper_bound) {
            sum_values += value;
            count++;
        }
    }

    // Return the mean of the filtered segment
    return count > 0 ? sum_values / count : 0;  // Ensure we don't divide by zero
}

void normalize(parlay::sequence<double> &signal) {
    size_t s_len = signal.size();

    double mean = std::accumulate(signal.begin(), signal.end(), 0.0) / s_len;
    double sq_sum = std::inner_product(signal.begin(), signal.end(), signal.begin(), 0.0);
    double std_dev = std::sqrt(sq_sum / s_len - mean * mean);

    size_t k = 0;
    for (size_t i = 0; i < s_len; ++i) {
        double norm_val = (signal[i] - mean) / std_dev;
        if (norm_val < 3 && norm_val > -3) {
            signal[k] = norm_val;
            k++;
        }
    }

    signal.resize(k);
}

parlay::sequence<double> generate_events(parlay::sequence<double> &signal, signal_segmenter_i &segmenter) {
    normalize(signal);
    auto peaks = segmenter.segment_signal(signal);
    uint32_t peak_size = peaks.size(), s_len = signal.size(), n_ev = 0;

    for (uint32_t pi = 0; pi < peak_size; ++pi)
        if (peaks[pi] > 0 && peaks[pi] < s_len) n_ev++;

    parlay::sequence<double> events(n_ev);

    uint32_t start_idx = 0, segment_length = 0;

    for (uint32_t pi = 0, i = 0; pi < peak_size && i < n_ev; pi++){
        if (!(peaks[pi] > 0 && peaks[pi] < s_len)) continue;

        segment_length = peaks[pi] - start_idx;
        events[i++] = calculate_mean_of_filtered_segment(signal.data() + start_idx, signal.data() + start_idx + segment_length);
        start_idx = peaks[pi];
    }

    return events;
}