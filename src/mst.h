#ifndef MST_H
#define MST_H

#include "data_type.h"
#include "graph.h"
#include "sparsehash/dense_hash_map"
#include "quad_tree.h"
#include "window_generator.h"
#include "hybrid.h"
#include "dynamic_mst.h"
#include "BS_thread_pool.hpp"
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <tuple>
#include <fstream>
#include <chrono>
#include <climits>

using google::dense_hash_map;

class MSTIndex {
public:
    TemporalGraph* graph_;
    bint time_window_sz;
    // indices
    std::vector<dense_hash_map<Edge, ExpireTimeRecord, edge_hash>> bl_idx_;
    std::vector<ExpireTimeRecord> expired_seq;
    std::vector<std::vector<Edge>> bl_keys_;
    std::vector<std::vector<ExpireTimeRecord>> bl_values_;
    std::vector<bint> op_t1_;
    std::vector<std::vector<bint>> op_t2_;
    std::vector<std::vector<std::vector<bint>>> op_s_;
    std::vector<std::vector<std::vector<std::vector<bint>>>> op_t_;
    std::vector<std::vector<std::vector<std::vector<std::vector<Edge>>>>> op_values_;
    std::unique_ptr<QuadTree> qt_idx_;
    std::unique_ptr<SegmentTree> hy_idx_;
    std::vector<Edge> cached_total_seq;
    std::vector<std::vector<int>> cached_basic_msts;

    explicit MSTIndex(TemporalGraph* graph) {
        graph_ = graph;
        time_window_sz = graph->get_num_timestamps();
    }
    static std::vector<Edge> Kruskal_MST(std::vector<Edge> edges, int num_nodes);
    static std::vector<Edge> Kruskal_MST_wo_sort(std::vector<Edge> edges, int num_nodes);
    static std::pair<std::vector<Edge>, std::vector<Edge>> special_Kruskal_MST(std::vector<Edge> edges, int num_nodes);
    void create_bl_idx();
    void merge(const std::vector<int>& nums1, const std::vector<int>& nums2, std::vector<int> &result);
    void create_bl_idx_plus();
    void create_bl_idx_pro();
    void update_bl_idx();
    void create_bl_idx_pro_parallel(int num_threads);
    void create_op_idx();
    void create_qt_idx();
    void create_hy_idx(bool xFirst);
    void qt_idx_count(std::string &index_file);
    std::vector<std::vector<Edge>> query_online(std::vector<TimeWindow> &time_windows) const;
    std::vector<std::vector<Edge>> query_online_global(std::vector<TimeWindow> &time_windows) const;
    std::vector<Edge> k_way_merge(std::vector<std::vector<Edge>>& sorted_lists);
    std::vector<std::vector<Edge>> query_online_k_way(std::vector<TimeWindow> &time_windows);
    std::vector<std::vector<Edge>> query_bl(std::vector<TimeWindow> &time_windows);
    std::vector<std::vector<Edge>> query_op(std::vector<TimeWindow> &time_windows);
    std::vector<std::vector<Edge>> query_qt(std::vector<TimeWindow> &time_windows);
    std::vector<std::vector<Edge>> query_hy(std::vector<TimeWindow> &time_windows);
    void validate(std::vector<TimeWindow> &time_windows, std::string index, std::string ratio);
    void save_bl_idx(std::string &index_file);
    void load_bl_idx(std::string &index_file);
    void save_op_idx(std::string &index_file);
    void load_op_idx(std::string &index_file);
    void save_qt_idx(std::string &index_file);
    void load_qt_idx(std::string &index_file);
    void save_hy_idx(std::string &index_file);
    void load_hy_idx(std::string &index_file);
    void delete_op_idx();
    void delete_qt_idx();
    void delete_hy_idx();
    template <typename T> void save_vector(std::ofstream &out, const std::vector<T> &vec);
    template <typename T> void save_nested_vector(std::ofstream &out, const std::vector<std::vector<T>> &nested_vec);
    template <typename T> void save_double_nested_vector(std::ofstream &out, const std::vector<std::vector<std::vector<T>>> &nested_vec);
    template <typename T> void save_triple_nested_vector(std::ofstream &out, const std::vector<std::vector<std::vector<std::vector<T>>>> &nested_vec);
    template <typename T> void save_quadruple_nested_vector(std::ofstream &out, const std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>> &nested_vec);
    template <typename T> void load_vector(std::ifstream &in, std::vector<T> &vec);
    template <typename T> void load_nested_vector(std::ifstream &in, std::vector<std::vector<T>> &nested_vec);
    template <typename T> void load_double_nested_vector(std::ifstream &in, std::vector<std::vector<std::vector<T>>> &nested_vec);
    template <typename T> void load_triple_nested_vector(std::ifstream &in, std::vector<std::vector<std::vector<std::vector<T>>>> &nested_vec);
    template <typename T> void load_quadruple_nested_vector(std::ifstream &in, std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>> &nested_vec);
};

#endif // MST_H