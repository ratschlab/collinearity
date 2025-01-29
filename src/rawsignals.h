//
// Created by Sayan Goswami on 23.01.2025.
//

#ifndef COLLINEARITY_RAWSIGNALS_H
#define COLLINEARITY_RAWSIGNALS_H

#include <utility>

#include "prelude.h"
#include "vector"

std::pair<int, std::vector<double>> load_pore_model(std::string &poremodel_file);

//parlay::sequence<float> sequence2squiggles(const std::string &sequence, int k, const std::vector<float> &levels);

template <typename T>
parlay::sequence<double> sequence2squiggles(const T &sequence, int k, const std::vector<double> &levels);

parlay::sequence<u1> quantize_signal(parlay::sequence<double> &signal, float diff=.35f, u1 quant_bit=4,
                                float fine_min=-2.0f, float fine_max=2.0f, float fine_range=.4f);

class signal_segmenter_i {
public:
    virtual ~signal_segmenter_i() {}
    virtual parlay::sequence<size_t> segment_signal(const parlay::sequence<double> &signal) = 0;
};

parlay::sequence<double> generate_events(parlay::sequence<double> &signal, signal_segmenter_i &segmenter);


struct ri_detect_t {
    const int DEF_PEAK_POS = -1;
    const float DEF_PEAK_VAL = __FLT_MAX__;
    const float THRESHOLD;
    const uint32_t WINDOW_LENGTH;

    parlay::sequence<double> tstat;
    uint32_t masked_to = 0;
    int peak_pos = -1;
    double peak_value = __FLT_MAX__;
    int valid_peak = 0;
    ri_detect_t(u4 window_length, float threshold): WINDOW_LENGTH(window_length), THRESHOLD(threshold) {}
    void reset() { masked_to = 0, peak_pos = DEF_PEAK_POS, peak_value = DEF_PEAK_VAL, valid_peak = 0; }
};

class tstat_segmenter_t : public signal_segmenter_i {
    std::vector<ri_detect_t> detectors;
public:
    tstat_segmenter_t() {
        detectors.emplace_back(3, 4.0);
        detectors.emplace_back(9, 3.5);
    }
    parlay::sequence<size_t> segment_signal(const parlay::sequence<double> &signal);
};

#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))


#endif //COLLINEARITY_RAWSIGNALS_H
