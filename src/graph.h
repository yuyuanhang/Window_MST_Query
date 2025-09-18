#ifndef GRAPH_H
#define GRAPH_H

#include "data_type.h"
#include <vector>
#include <unordered_set>
#include <map>

class TemporalGraph {
public:
    int nb_tbase;
    int nb_tstream;
    explicit TemporalGraph(const std::vector<Edge>& temporal_edges, float edge_ratio = 1);

    void group_edges_by_timestamp();
    std::vector<Edge> get_temporal_edges(bint timestamp) { 
        if (timestamp < num_timestamps_)
            return edges_by_timestamp_[timestamp];
        else 
            return std::vector<Edge>{};
    }
    [[nodiscard]] int get_num_nodes() const { return num_nodes_; }
    [[nodiscard]] bint get_num_timestamps() const { return num_timestamps_; }
    [[nodiscard]] bint get_num_edges() const { return num_edges_; }
    [[maybe_unused]] std::vector<Edge> edges() { return edges_; }
    [[maybe_unused]] void print_info() const;

private:
    std::vector<Edge> edges_;
    std::vector<std::vector<Edge>> edges_by_timestamp_;

    int num_nodes_;
    bint num_edges_;
    bint num_timestamps_;
};

#endif // GRAPH_H