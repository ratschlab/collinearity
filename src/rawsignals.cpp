//
// Created by Sayan Goswami on 23.01.2025.
//

#include "collinearity.h"

u1 dynamic_quantize(double signal, float fine_min, float fine_max, float fine_range, u4 n_buckets) {
    // Total range for normalization
    float minVal = -3.0, maxVal = 3.0;
    float range = maxVal - minVal;
    float coarse_coef1 = (1-fine_range)/2;
    float coarse_coef2 = fine_range + coarse_coef1;

    // Normalize the signal to [0, 1]
    double normalized = (signal - minVal) / range;

    float a = (fine_min - minVal) / range;
    float b = (fine_max - minVal) / range;

    // Conditional quantization based on the segment
    double quantized = fine_max;
    if (signal >= fine_min && signal <= fine_max) {
        // Within [fine_min, fine_max], map to a sub-range [a, b] in [0, 1],
        //then scale to [0, fine_range] for finer granularity
        quantized = fine_range * ((normalized - a) / (b - a));
    }
    else {
        // Outside [fine_min, fine_max], split the rest of [0, 1] into two and map accordingly
        if (normalized < 0.5) quantized = fine_range + coarse_coef1 * normalized; // Coarser granularity
        else quantized = coarse_coef2 + coarse_coef1 * normalized; // Coarser granularity
    }

    // Map the quantized value back to the range [0, 2^n_buckets - 1]
    u1 quantizedValue = (u1)(quantized * (n_buckets-1));

    return quantizedValue;
}

parlay::sequence<u1> quantize_signal(parlay::sequence<double> &signal, float diff, u1 quant_bit,
                                float fine_min, float fine_max, float fine_range) {
    const auto slen = signal.size();
    parlay::sequence<u1> quantized(slen);

    int n_buckets = 1 << quant_bit;
    int mask_quant_bit = (1 << quant_bit) - 1;

    int f_pos = 0;
    int l_sigpos = 0;  // last signal position
    int nqval = 0; // numbver of quantized values

    // First quantization is done here
    l_sigpos = f_pos;
    quantized[nqval++] = dynamic_quantize(signal[f_pos], fine_min, fine_max, fine_range, n_buckets) & mask_quant_bit;

    for (f_pos = 1; f_pos < slen; f_pos++) {
        if (std::abs(signal[f_pos] - signal[l_sigpos]) < diff) {
            continue;
        }

        l_sigpos = f_pos;
        quantized[nqval++] = dynamic_quantize(signal[f_pos], fine_min, fine_max, fine_range, n_buckets) & mask_quant_bit;
    }

    quantized.resize(nqval);
    return quantized;
}