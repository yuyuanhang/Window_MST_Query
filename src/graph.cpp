#include "graph.h"
#include <iostream>
#include <algorithm>
#include <random>

TemporalGraph::TemporalGraph(const std::vector<Edge>& temporal_edges, float edge_ratio) {
    edges_ = temporal_edges;

    std::sort(edges_.begin(), edges_.end(), [](const Edge& a, const Edge& b) {
        return a.timestamp < b.timestamp;
    });

    size_t original_size = edges_.size();
    auto num_elements_to_copy = static_cast<size_t>(std::ceil(static_cast<double>(original_size) * edge_ratio));
    edges_ = std::vector<Edge>(edges_.begin(), edges_.begin() + num_elements_to_copy);
    num_edges_ = static_cast<bint>(edges_.size());

    num_nodes_ = 0;
    num_timestamps_ = 0;
    for (auto edge : edges_) {
        if (edge.src_id >= num_nodes_) {
            num_nodes_ = edge.src_id + 1;
        }
        if (edge.dest_id >= num_nodes_) {
            num_nodes_ = edge.dest_id + 1;
        }
        if (edge.timestamp >= num_timestamps_) { num_timestamps_ = edge.timestamp+1; }
    }

    std::cout << "|V|: " << num_nodes_ << std::endl;
    std::cout << "|E|: " << num_edges_ << std::endl;
    std::cout << "t_max: " << num_timestamps_ << std::endl;
}

void TemporalGraph::group_edges_by_timestamp() {
    edges_by_timestamp_.resize(num_timestamps_, std::vector<Edge>());
    for (auto &edge : edges_) {
        edges_by_timestamp_[edge.timestamp].emplace_back(edge);
    }
}

[[maybe_unused]] void TemporalGraph::print_info() const {
    std::cout << "|V|: " << num_nodes_ << std::endl;
    std::cout << "|E|: " << num_edges_ << std::endl;
    for (int i = 0; i < 5; i++) {
        std::cout << "Edge " << i << ": " << 
            edges_[i].src_id << " " << edges_[i].dest_id << " " <<
            edges_[i].weight << " " << edges_[i].timestamp << std::endl;
    }
    for (int i = num_edges_ - 5; i <= num_edges_ - 1; i++) {
        std::cout << "Edge " << i << ": " << 
            edges_[i].src_id << " " << edges_[i].dest_id << " " <<
            edges_[i].weight << " " << edges_[i].timestamp << std::endl;
    }
    std::cout << "t_max: " << num_timestamps_ << std::endl;
}