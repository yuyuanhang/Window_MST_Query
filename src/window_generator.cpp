//
// Created by Yuanhang Yu on 20/11/2023.
//

#include "window_generator.h"

std::vector<TimeWindow> TimeWindowGenerator::generate(int count, float window_ratio) {
    std::cout << "generate time windows..." << std::endl;
    std::vector<TimeWindow> windows;

    int window_size = std::floor(static_cast<float>(num_timestamp_-1) * window_ratio);
    for (int i = 0; i < count; ++i) {
        int start_time = rand_range(0, num_timestamp_-(window_size+1));

        TimeWindow window;
        window.first = start_time;
        window.second = start_time + window_size;
        windows.push_back(window);
    }

    windows_ = windows;
    return windows;
}

int TimeWindowGenerator::rand_range(int min, int max) {
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator_);
};

void TimeWindowGenerator::save_windows() {
    std::cout << "save time windows..." << std::endl;
    std::ofstream file(window_file_, std::ios::binary);

    if (!file) {
        throw std::runtime_error("can not open given file");
    }

    size_t size = windows_.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));

    for (const auto& window : windows_) {
        file.write(reinterpret_cast<const char*>(&window.first), sizeof(window.first));
        file.write(reinterpret_cast<const char*>(&window.second), sizeof(window.second));
    }

    file.close();
}

void TimeWindowGenerator::load_windows() {
    std::cout << "load time windows..." << std::endl;
    std::ifstream file(window_file_, std::ios::binary);

    if (!file) {
        throw std::runtime_error("can not open window file");
    }

    size_t size;
    file.read(reinterpret_cast<char*>(&size), sizeof(size));

    windows_.clear();
    windows_.resize(size);
    for (auto& window : windows_) {
        file.read(reinterpret_cast<char*>(&window.first), sizeof(window.first));
        file.read(reinterpret_cast<char*>(&window.second), sizeof(window.second));
    }

    file.close();
}

[[maybe_unused]] void TimeWindowGenerator::print_windows() {
    for (int i = 0; i < 5; i++) {
        std::cout << "Time window " << i << ": [" <<
                  windows_[i].first << ", " << windows_[i].second << "]" << std::endl;
    }

    for (auto i = windows_.size() - 5; i <= windows_.size() - 1; i++) {
        std::cout << "Time window " << i << ": [" <<
                  windows_[i].first << ", " << windows_[i].second << "]" << std::endl;
    }
}