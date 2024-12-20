//
// Created by Yuanhang Yu on 20/11/2023.
//

#ifndef WINDOW_GENERATOR_H
#define WINDOW_GENERATOR_H
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <set>
#include "data_type.h"

class TimeWindowGenerator {
public:
    int num_timestamp_;
    std::default_random_engine generator_;
    std::vector<TimeWindow> windows_;
    std::string window_file_;

    TimeWindowGenerator(int num_timestamp, std::string &window_file) : num_timestamp_(num_timestamp), window_file_(window_file) {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator_.seed(seed);
    }

    std::vector<TimeWindow> generate(int count, float window_ratio);
    int rand_range(int min, int max);
    void save_windows();
    void load_windows();

    [[maybe_unused]] void print_windows();
};


#endif //WINDOW_GENERATOR_H
