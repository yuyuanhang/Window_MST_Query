#include "mst.h"
#include <iostream>
#include <algorithm>
#if defined(_WIN32) || defined(_WIN64)
#define OS_WINDOWS
#elif defined(__APPLE__) || defined(__MACH__)
#define OS_MAC
#elif defined(__linux__)
#define OS_LINUX
#else
#error "Unknown operating system"
#endif
#ifdef OS_WINDOWS
#include <windows.h>
#include <psapi.h>

void print_peak_memory_usage() {
    HANDLE hProcess = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS pmc;

    if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) {
        std::cout << "Peak memory usage: " << pmc.PeakWorkingSetSize / 1024 / 1024 << " MB" << std::endl;
    } else {
        std::cerr << "Failed to get memory info" << std::endl;
    }
}
#endif

#ifdef OS_MAC
#include <mach/mach.h>

void print_peak_memory_usage() {
    mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
        std::cout << "Peak memory usage: " << info.resident_size_max / 1024 / 1024 << " MB" << std::endl;
    } else {
        std::cerr << "Failed to get memory info" << std::endl;
    }
}
#endif

#ifdef OS_LINUX
#include <sys/resource.h>

void print_peak_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "Peak memory usage: " << usage.ru_maxrss / 1024 << " MB" << std::endl;
}
#endif

std::vector<Edge> MSTIndex::Kruskal_MST(std::vector<Edge> edges, int num_nodes) {
    std::vector<Edge> result;

    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        if (a.weight != b.weight) {
            return a.weight < b.weight;
        } else if (a.timestamp != b.timestamp) {
            return a.timestamp > b.timestamp;
        } else if (a.src_id != b.src_id) {
            return a.src_id < b.src_id;
        } else {
            return a.dest_id < b.dest_id;
        }
    });

    UnionFind uf(num_nodes);
    uf.clear(num_nodes);

    for (auto edge : edges) {
        int u = edge.src_id;
        int v = edge.dest_id;
        int rootX, rootY;

        if (!uf.connected(u, v, rootX, rootY)) {
            result.push_back(edge);
            uf.unite(rootX, rootY);
        }
    }

    return result;
}

std::vector<Edge> MSTIndex::Kruskal_MST_wo_sort(std::vector<Edge> edges, int num_nodes)
{
    std::vector<Edge> result;

    UnionFind uf(num_nodes);
    uf.clear(num_nodes);

    for (auto edge : edges) {
        int u = edge.src_id;
        int v = edge.dest_id;
        int rootX, rootY;

        if (!uf.connected(u, v, rootX, rootY)) {
            result.push_back(edge);
            uf.unite(rootX, rootY);
        }
    }

    return result;
}

std::pair<std::vector<Edge>, std::vector<Edge>> MSTIndex::special_Kruskal_MST(std::vector<Edge> edges, int num_nodes) {
    std::vector<Edge> result;
    std::vector<Edge> invalid;

    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        if (a.weight != b.weight) {
            return a.weight < b.weight;
        } else if (a.timestamp != b.timestamp) {
            return a.timestamp > b.timestamp;
        } else if (a.src_id != b.src_id) {
            return a.src_id < b.src_id;
        } else {
            return a.dest_id < b.dest_id;
        }
    });

    UnionFind uf(num_nodes);
    uf.clear(num_nodes);

    for (auto edge : edges) {
        int u = edge.src_id;
        int v = edge.dest_id;
        int rootX, rootY;

        if (!uf.connected(u, v, rootX, rootY)) {
            result.push_back(edge);
            uf.unite(rootX, rootY);
        } else {
            invalid.push_back(edge);
        }
    }

    return std::make_pair(result, invalid);
}

void MSTIndex::create_bl_idx() {
    auto start = std::chrono::high_resolution_clock::now();
    bl_idx_.clear();
    bl_idx_.resize(time_window_sz);
    std::vector<std::vector<Edge>> basic_msts;
    basic_msts.resize(time_window_sz);

    auto basic_edges = graph_->get_temporal_edges(0);
    auto mst = Kruskal_MST(basic_edges, graph_->get_num_nodes());
    basic_msts[0] = mst;
    auto seg = dense_hash_map<Edge, ExpireTimeRecord, edge_hash>();
    seg.set_empty_key(Edge(0, 0, 0, 0));
    for (auto e : mst) {
        seg[e] = ExpireTimeRecord();
        seg[e].emplace_back(std::make_pair(-1, 0));
    }
    bl_idx_[0] = seg;

    for (auto i = 1; i < time_window_sz; i++) {
        basic_edges = graph_->get_temporal_edges(i);
        mst = Kruskal_MST(basic_edges, graph_->get_num_nodes());
        basic_msts[i] = mst;
        seg = dense_hash_map<Edge, ExpireTimeRecord, edge_hash>();
        seg.set_empty_key(Edge(0, 0, 0, 0));
        for (auto e : mst) {
            seg[e] = ExpireTimeRecord();
        }
        bl_idx_[i] = seg;

        size_t cnt_expired = 0;
        auto target_expired = basic_msts[i].size();
        for (auto j = i - 1; j >= 0; j--) {
            auto added_mst = basic_msts[j];
            mst.insert(mst.end(), added_mst.begin(), added_mst.end());
            auto p = special_Kruskal_MST(mst, graph_->get_num_nodes());
            for (auto e : p.second) {
                auto record = bl_idx_[e.timestamp][e];
                if (record.empty() || record.back().first != j) {
                    bl_idx_[e.timestamp][e].emplace_back(std::make_pair(j, i));
                }
                if (e.timestamp == i) {
                    cnt_expired++;
                }
            }

            mst = p.first;
            if (j == 0) {
                for (auto e : mst) {
                    auto record = bl_idx_[e.timestamp][e];
                    if (record.empty() || record.back().first != -1) {
                        bl_idx_[e.timestamp][e].emplace_back(std::make_pair(-1, i));
                    }
                }
            }
            if (cnt_expired == target_expired) {
                break;
            }
        }
    }

    bl_keys_.resize(time_window_sz);
    bl_values_.resize(time_window_sz);
    for (auto i = 0; i < time_window_sz; i++) {
        for (auto &p : bl_idx_[i]) {
            bl_keys_[i].push_back(p.first);
            bl_values_[i].push_back(p.second);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "bl index construction time: " << (double)duration.count()/1000 << " ms" << std::endl;
    print_peak_memory_usage();
}

void MSTIndex::merge(const std::vector<int>& nums1, const std::vector<int>& nums2, std::vector<int> &result) {
    auto size_a = nums1.size();
    auto size_b = nums2.size();
    result.reserve(size_a+size_b);
    int i = 0, j = 0;

    while (i < size_a && j < size_b) {
        auto ele_a = nums1[i];
        auto ele_b = nums2[j];
        if (ele_a < ele_b) {
            result.push_back(ele_a);
            i++;
        } else {
            result.push_back(ele_b);
            j++;
        }
    }
    if (i < size_a) {
        result.insert(result.end(), nums1.begin() + i, nums1.end());
    }
    if (j < size_b) {
        result.insert(result.end(), nums2.begin() + j, nums2.end());
    }
}

void MSTIndex::create_bl_idx_plus() {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Edge> total_seq; // edges themselves

    for (auto i = 0; i < time_window_sz; i++) {
        // compute basic mst $T_{[i, i]}$
        auto basic_edges = graph_->get_temporal_edges(i);
        auto mst = Kruskal_MST(basic_edges, graph_->get_num_nodes());

        // store edges
        total_seq.insert(total_seq.end(), mst.begin(), mst.end());
    }

    // sort total_seq for subsequent mst computation
    std::sort(total_seq.begin(), total_seq.end(), [](const Edge& a, const Edge& b) {
        if (a.weight != b.weight) {
            return a.weight < b.weight;
        } else if (a.timestamp != b.timestamp) {
            return a.timestamp > b.timestamp;
        } else if (a.src_id != b.src_id) {
            return a.src_id < b.src_id;
        } else {
            return a.dest_id < b.dest_id;
        }
    });

    // an inverted index
    std::vector<std::vector<int>> ordered_edges; // basic msts
    ordered_edges.resize(time_window_sz);
    int cnt = 0;
    for (auto &edge : total_seq) {
        ordered_edges[edge.timestamp].push_back(cnt);
        cnt++;
    }

    // ExpireTimeRecord vector
    expired_seq.resize(total_seq.size());
    for (auto idx : ordered_edges[0]) {
        expired_seq[idx].emplace_back(-1, 0);
    }

    // Union
    auto num_nodes = graph_->get_num_nodes();
    UnionFind uf(num_nodes);
    for (auto i = 1; i < time_window_sz; i++) {
        auto mst_size = ordered_edges[i].size(); // for early termination
        std::vector<int> local_edges; // $T_{[i, i]}$, and becomes $T_{[j, i]}$
        local_edges = ordered_edges[i];

        size_t cnt_expired = 0; // counter for expired edges
        for (auto j = i - 1; j >= 0; j--) {
            std::vector<int> added_edges(ordered_edges[j]); // $T_{[j, j]}$

            // merge two sorted lists to get $T_{[j, i]}$
            std::vector<int> expired_edges; // for $T_{[j, j]}$ when j != 0
            std::vector<int> merged_edges; // for $T_{[j, j]}$ when j == 0
            merge(local_edges, added_edges, merged_edges);

            local_edges.clear();
            uf.clear(num_nodes);
            for (auto k = 0; k < merged_edges.size(); k++) {
                auto idx = merged_edges[k];
                auto crt_e = total_seq[idx];
                int u = crt_e.src_id;
                int v = crt_e.dest_id;
                int rootX = 0, rootY = 0;

                if (!uf.connected(u, v, rootX, rootY)) {
                    uf.unite(rootX, rootY);
                    local_edges.push_back(idx);
                } else {
                    expired_edges.push_back(idx);
                }
            }

            for (auto& idx : expired_edges) {
                auto record = expired_seq[idx];
                if (record.empty() || record.back().first < j) {
                    expired_seq[idx].emplace_back(j, i);
                }
                if (total_seq[idx].timestamp == i) {
                    cnt_expired++;
                }
            }

            if (j == 0) {
                for (auto& idx : local_edges) {
                    auto record = expired_seq[idx];
                    if (record.empty()) {
                        expired_seq[idx].emplace_back(-1, i);
                    }
                }
            }
            if (cnt_expired == mst_size) {
                break;
            }
        }
    }

    bl_keys_.resize(time_window_sz);
    bl_values_.resize(time_window_sz);
    cnt = 0;
    for (auto& edge_list : ordered_edges) {
        for (auto &idx : edge_list) {
            bl_keys_[cnt].push_back(total_seq[idx]);
            bl_values_[cnt].push_back(expired_seq[idx]);
        }
        cnt++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "bl index construction time (plus): " << (double)duration.count()/1000 << " ms" << std::endl;
    print_peak_memory_usage();
}

std::vector<Edge> merge_sorted_edges(const std::vector<Edge>& arr1, const std::vector<Edge>& arr2) {
    std::vector<Edge> merged;
    size_t i = 0, j = 0;
    size_t n1 = arr1.size();
    size_t n2 = arr2.size();

    /*while (i < n1 && j < n2) {
        if (edge_order({arr1[i], 0}, {arr2[j], 0})) {
            merged.push_back(arr1[i]);
            i++;
        } else {
            merged.push_back(arr2[j]);
            j++;
        }
    }*/

    // Add remaining elements from arr1
    while (i < n1) {
        merged.push_back(arr1[i]);
        i++;
    }

    // Add remaining elements from arr2
    while (j < n2) {
        merged.push_back(arr2[j]);
        j++;
    }

    return merged;
}

void MSTIndex::create_bl_idx_pro() {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Edge> total_seq;
    std::vector<int> threshold;
    for (auto i = 0; i < time_window_sz; i++) {
        auto basic_edges = graph_->get_temporal_edges(i);
        auto mst = Kruskal_MST(basic_edges, graph_->get_num_nodes());
        threshold.push_back(mst.size());
        if (i == 0) {
            total_seq = mst;
        } else {
            total_seq = merge_sorted_edges(total_seq, mst);
        }
    }

    std::vector<std::vector<int>> basic_msts;
    basic_msts.resize(time_window_sz);
    int cnt = 0;
    for (auto &edge : total_seq) {
        basic_msts[edge.timestamp].push_back(cnt);
        cnt++;
    }

    // ExpireTimeRecord vector
    expired_seq.resize(total_seq.size());
    for (auto idx : basic_msts[0]) {
        expired_seq[idx].emplace_back(-1, 0);
    }

    // Union
    auto num_nodes = graph_->get_num_nodes();
    for (auto i = 1; i < time_window_sz; i++) {
        DMST::DynamicMST dMST(num_nodes);
        dMST.init();
        std::unordered_map<DMST::Edge, int, DMST::EdgeHash> edgeMap;

        int cnt_expired = 0;
        for (auto &e_idx : basic_msts[i]) {
            Edge edge = total_seq[e_idx];
            DMST::Edge removed_e;
            dMST.insert_edge(edge.src_id, edge.dest_id, edge.weight, removed_e);
            DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
            edgeMap[added_edge] = e_idx;
        }

        for (auto j = i - 1; j >= 0; j--) {
            for (auto &e_idx : basic_msts[j]) {
                Edge edge = total_seq[e_idx];
                DMST::Edge removed_e;
                
                if (expired_seq[e_idx].back().first < j) {
                    auto res = dMST.insert_edge(edge.src_id, edge.dest_id, edge.weight, removed_e);
                    if (res == 1) { // fail
                        expired_seq[e_idx].emplace_back(j, i);
                        if (total_seq[e_idx].timestamp == i) {
                            cnt_expired++;
                        }
                    } else if (res == 2) {
                        if (removed_e.u > removed_e.v) { std::swap(removed_e.u, removed_e.v); }
                        auto removed_e_idx = edgeMap[removed_e];
                        DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
                        edgeMap[added_edge] = e_idx;
                        if (expired_seq[removed_e_idx].empty()) {
                            expired_seq[removed_e_idx].emplace_back(j, i);
                        }
                        else if (expired_seq[removed_e_idx].back().first < j) {
                            expired_seq[removed_e_idx].emplace_back(j, i);
                        }
                        edgeMap.erase(removed_e);
                        if (total_seq[removed_e_idx].timestamp == i) {
                            cnt_expired++;
                        }
                    } else {
                        DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
                        edgeMap[added_edge] = e_idx;
                    }
                }
            }

            if (j == 0) {
                for (auto& e_idx : basic_msts[i]) {
                    if (expired_seq[e_idx].empty()) {
                        expired_seq[e_idx].emplace_back(-1, i);
                    }
                }
            }
            if (cnt_expired == threshold[i]) {
                break;
            }
        }
    }

    bl_keys_.resize(time_window_sz);
    bl_values_.resize(time_window_sz);
    cnt = 0;
    for (auto& edge_list : basic_msts) {
        for (auto &idx : edge_list) {
            bl_keys_[cnt].push_back(total_seq[idx]);
            bl_values_[cnt].push_back(expired_seq[idx]);
        }
        cnt++;
    }
    cached_total_seq = std::move(total_seq);
    cached_basic_msts = std::move(basic_msts);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "bl index construction time (pro): " << (double)duration.count()/1000 << " ms" << std::endl;
    print_peak_memory_usage();
}

void MSTIndex::update_bl_idx()
{
    auto basic_edges = graph_->get_temporal_edges(time_window_sz);
    auto mst = Kruskal_MST(basic_edges, graph_->get_num_nodes());
    int threshold = static_cast<int>(mst.size());
    cached_basic_msts.resize(time_window_sz + 1);
    for (int i = 0; i < mst.size(); i++)
    {
        cached_basic_msts[time_window_sz].push_back(i + cached_total_seq.size());
    }

    cached_total_seq = merge_sorted_edges(cached_total_seq, mst);

    // ExpireTimeRecord vector
    expired_seq.resize(cached_total_seq.size());

    // Union
    auto num_nodes = graph_->get_num_nodes();
    DMST::DynamicMST dMST(num_nodes);
    dMST.init();
    std::unordered_map<DMST::Edge, int, DMST::EdgeHash> edgeMap;

    int cnt_expired = 0;
    for (auto &e_idx : cached_basic_msts[time_window_sz]) {
        Edge edge = cached_total_seq[e_idx];
        DMST::Edge removed_e;
        dMST.insert_edge(edge.src_id, edge.dest_id, edge.weight, removed_e);
        DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
        edgeMap[added_edge] = e_idx;
    }

    for (auto i = time_window_sz - 1; i >= 0; i--) {
        for (auto &e_idx : cached_basic_msts[i]) {
            Edge edge = cached_total_seq[e_idx];
            DMST::Edge removed_e;

            if (expired_seq[e_idx].back().first < i) {
                auto res = dMST.insert_edge(edge.src_id, edge.dest_id, edge.weight, removed_e);
                if (res == 1) { // fail
                    expired_seq[e_idx].emplace_back(i, time_window_sz);
                    if (cached_total_seq[e_idx].timestamp == time_window_sz) {
                        cnt_expired++;
                    }
                } else if (res == 2) {
                    if (removed_e.u > removed_e.v) { std::swap(removed_e.u, removed_e.v); }
                    auto removed_e_idx = edgeMap[removed_e];
                    DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
                    edgeMap[added_edge] = e_idx;
                    if (expired_seq[removed_e_idx].empty()) {
                        expired_seq[removed_e_idx].emplace_back(i, time_window_sz);
                    }
                    else if (expired_seq[removed_e_idx].back().first < i) {
                        expired_seq[removed_e_idx].emplace_back(i, time_window_sz);
                    }
                    edgeMap.erase(removed_e);
                    if (cached_total_seq[removed_e_idx].timestamp == time_window_sz) {
                        cnt_expired++;
                    }
                } else {
                    DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
                    edgeMap[added_edge] = e_idx;
                }
            }
        }

        if (i == 0) {
            for (auto& e_idx : cached_basic_msts[time_window_sz]) {
                if (expired_seq[e_idx].empty()) {
                    expired_seq[e_idx].emplace_back(-1, i);
                }
            }
        }

        if (cnt_expired == threshold) {
            break;
        }
    }

    bl_keys_.resize(time_window_sz + 1);
    bl_values_.resize(time_window_sz + 1);
    for (auto i = 0; i < time_window_sz; i++) {
        int cnt = 0;
        for (auto &idx : cached_basic_msts[i]) {
            bl_values_[i][cnt] = expired_seq[idx];
            cnt++;
        }
    }
    for (auto &idx : cached_basic_msts[time_window_sz])
    {
        bl_keys_[time_window_sz].push_back(cached_total_seq[idx]);
        bl_values_[time_window_sz].push_back(expired_seq[idx]);
    }
    time_window_sz++;
}

void MSTIndex::create_bl_idx_pro_parallel(int num_threads) {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Edge> total_seq;
    std::vector<int> threshold;
    std::vector<std::vector<Edge>> mst_cache(time_window_sz);
    threshold.resize(time_window_sz);

    BS::thread_pool pool(num_threads);

    const BS::multi_future<void> loop_future = pool.submit_loop(0, time_window_sz,
        [&](int i) {
            auto basic_edges = graph_->get_temporal_edges(i);
            mst_cache[i] = Kruskal_MST(basic_edges, graph_->get_num_nodes());
            threshold[i] = static_cast<int>(mst_cache[i].size());
        }
    );

    loop_future.wait();

    for (int i = 0; i < time_window_sz; i++) {
        if (i == 0) {
            total_seq = std::move(mst_cache[i]);
        } else {
            total_seq.insert(total_seq.end(), mst_cache[i].begin(), mst_cache[i].end());
        }
    }
    mst_cache.clear();

    std::vector<std::vector<int>> basic_msts(time_window_sz);
    {
        int idx = 0;
        for (auto &e : total_seq) {
            basic_msts[e.timestamp].push_back(idx++);
        }
    }

    expired_seq.clear();
    expired_seq.resize(total_seq.size());
    for (int idx : basic_msts[0]) {
        expired_seq[idx].emplace_back(-1, 0);
    }

    auto num_nodes = graph_->get_num_nodes();

    int batch_size = num_threads;

    for (int i = 1; i < time_window_sz; i += batch_size) {
        int task_begin = i;
        int task_end = min(time_window_sz, i + batch_size);

        std::vector<std::unordered_map<int, std::pair<bint, bint>>> record_cache(task_end - task_begin);
        BS::multi_future<void> loop_compute;
        loop_compute.reserve(task_end - task_begin);
        for (int j = task_end - 1; j >= task_begin; j--)
        {
            loop_compute.push_back(pool.submit_task([&, j]
            {
                int local_id = (j - 1) % batch_size;

                DMST::DynamicMST dMST(num_nodes);
                dMST.init();
                std::unordered_map<DMST::Edge, int, DMST::EdgeHash> edgeMap;

                int cnt_expired = 0;
                for (auto &e_idx : basic_msts[j]) {
                    Edge edge = total_seq[e_idx];
                    DMST::Edge removed_e;
                    dMST.insert_edge(edge.src_id, edge.dest_id, edge.weight, removed_e);
                    DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
                    edgeMap[added_edge] = e_idx;
                }

                for (auto k = j - 1; k >= 0; k--) {
                    for (auto &e_idx : basic_msts[k]) {
                        Edge edge = total_seq[e_idx];
                        DMST::Edge removed_e;

                        if (expired_seq[e_idx].empty() || expired_seq[e_idx].back().first < k) {
                            auto res = dMST.insert_edge(edge.src_id, edge.dest_id, edge.weight, removed_e);
                            if (res == 1) { // fail
                                record_cache[local_id][e_idx] = {k, j};

                                if (edge.timestamp == j) {
                                    cnt_expired++;
                                }
                            } else if (res == 2) {
                                if (removed_e.u > removed_e.v) { std::swap(removed_e.u, removed_e.v); }
                                auto removed_e_idx = edgeMap[removed_e];

                                record_cache[local_id][removed_e_idx] = {k, j};
                                edgeMap.erase(removed_e);

                                DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
                                edgeMap[added_edge] = e_idx;

                                if (total_seq[removed_e_idx].timestamp == j) {
                                    cnt_expired++;
                                }
                            } else {
                                DMST::Edge added_edge(edge.src_id, edge.dest_id, edge.weight);
                                edgeMap[added_edge] = e_idx;
                            }
                        }
                    }

                    if (k == 0) {
                        for (auto &e_idx : basic_msts[j]) {
                            if (!record_cache[local_id].count(e_idx)) {
                                record_cache[local_id][e_idx] = {-1, j};
                            }
                        }
                    }
                    if (cnt_expired == threshold[j]) break;
                }
            }));
        }

        loop_compute.wait();

        std::vector<std::vector<std::unordered_map<int, std::pair<bint, bint>>>> cache_split(
            batch_size, std::vector<std::unordered_map<int, std::pair<bint, bint>>>(task_end - task_begin)
        );

        for (int j = 0; j < task_end - task_begin; j++) {
            for (const auto& [key, val] : record_cache[j]) {
                int tid = key % batch_size;
                cache_split[tid][j][key] = val;
            }
        }

        BS::multi_future<void> loop_update;
        loop_update.reserve(batch_size);

        int actual_batch_size = task_end - task_begin;
        for (int j = 0; j < batch_size; j++)
        {
            loop_update.push_back(pool.submit_task([&, j, actual_batch_size]
            {
                for (int k = 0; k < actual_batch_size; k++)
                {
                    for (auto &p : cache_split[j][k])
                    {
                        int e_idx = p.first;
                        if (expired_seq[e_idx].empty())
                        {
                            expired_seq[e_idx].emplace_back(p.second);
                        } else
                        {
                            if (p.second.first < expired_seq[e_idx].back().first)
                            {
                                cerr << "violate lemma" << endl;
                                exit(0);
                            } else if (p.second.first > expired_seq[e_idx].back().first)
                            {
                                expired_seq[e_idx].emplace_back(p.second);
                            }
                        }
                    }
                }
            }));
        }

        loop_update.wait();
        cache_split.clear();
        record_cache.clear();
    }

    bl_keys_.assign(time_window_sz, {});
    bl_values_.assign(time_window_sz, {});

    BS::multi_future<void> write;
    for (int t = 0; t < time_window_sz; ++t) {
        write.push_back(pool.submit_task([&, t]
        {
            for (int idx : basic_msts[t]) {
                bl_keys_[t].push_back(total_seq[idx]);
                bl_values_[t].push_back(expired_seq[idx]);
            }
        }));
    }
    write.wait();

    auto end = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "bl index construction time (pro && parallel): "
              << (double)dur.count() / 1000 << " ms" << std::endl;
    print_peak_memory_usage();
}

void MSTIndex::create_op_idx() {
    auto start = std::chrono::high_resolution_clock::now();
    //create_bl_idx();
    //create_bl_idx_plus();
    // (s, t1, t2, t)
    std::vector<Quadruple> items;
    for (auto i = 0; i < time_window_sz; i++) {
        for (auto j = 0; j < bl_keys_[i].size(); j++) {
            for (auto k = 0; k < bl_values_[i][j].size(); k++) {
                if (k == bl_values_[i][j].size() - 1) {
                    items.emplace_back(bl_values_[i][j][k].first, bl_values_[i][j][k].second, INT_MAX, bl_keys_[i][j]);
                } else {
                    items.emplace_back(bl_values_[i][j][k].first, bl_values_[i][j][k].second, bl_values_[i][j][k+1].second, bl_keys_[i][j]);
                }
            }
        }
    }
    std::sort(items.begin(), items.end(), compareQuadruples);

    // initialize
    bint target_t1 = std::get<1>(items[0]);
    bint target_t2 = std::get<2>(items[0]);
    bint target_s = std::get<0>(items[0]);
    Edge e = std::get<3>(items[0]);
    bint target_t = e.timestamp;
    size_t t1_size = 0;
    size_t t2_size;
    size_t s_size;
    size_t t_size;
    op_t1_.push_back(target_t1);
    t1_size++;
    op_t2_.resize(t1_size);
    op_t2_[t1_size-1].push_back(target_t2);
    t2_size = 1;
    op_s_.resize(t1_size);
    op_s_[t1_size-1].resize(t2_size);
    op_s_[t1_size-1][t2_size-1].push_back(target_s);
    s_size = 1;
    op_t_.resize(t1_size);
    op_t_[t1_size-1].resize(t2_size);
    op_t_[t1_size-1][t2_size-1].resize(s_size);
    op_t_[t1_size-1][t2_size-1][s_size-1].push_back(target_t);
    t_size = 1;
    op_values_.resize(t1_size);
    op_values_[t1_size-1].resize(t2_size);
    op_values_[t1_size-1][t2_size-1].resize(s_size);
    op_values_[t1_size-1][t2_size-1][s_size-1].resize(t_size);
    op_values_[t1_size-1][t2_size-1][s_size-1][t_size-1].push_back(e);

    for (auto i = 1; i < items.size(); i++) {
        bint t1 = std::get<1>(items[i]);
        bint t2 = std::get<2>(items[i]);
        bint s = std::get<0>(items[i]);
        e = std::get<3>(items[i]);
        bint t = e.timestamp;

        if (t1 != target_t1) {
            op_t1_.push_back(t1);
            t1_size++;
            op_t2_.resize(t1_size);
            op_t2_[t1_size-1].push_back(t2);
            t2_size = 1;
            op_s_.resize(t1_size);
            op_s_[t1_size-1].resize(t2_size);
            op_s_[t1_size-1][t2_size-1].push_back(s);
            s_size = 1;
            op_t_.resize(t1_size);
            op_t_[t1_size-1].resize(t2_size);
            op_t_[t1_size-1][t2_size-1].resize(s_size);
            op_t_[t1_size-1][t2_size-1][s_size-1].push_back(t);
            t_size = 1;
            op_values_.resize(t1_size);
            op_values_[t1_size-1].resize(t2_size);
            op_values_[t1_size-1][t2_size-1].resize(s_size);
            op_values_[t1_size-1][t2_size-1][s_size-1].resize(t_size);
            op_values_[t1_size-1][t2_size-1][s_size-1][t_size-1].push_back(e);
            target_t1 = t1;
            target_t2 = t2;
            target_s = s;
            target_t = t;
            continue;
        }

        if (t2 != target_t2) {
            op_t2_[t1_size-1].push_back(t2);
            t2_size++;
            op_s_[t1_size-1].resize(t2_size);
            op_s_[t1_size-1][t2_size-1].push_back(s);
            s_size = 1;
            op_t_[t1_size-1].resize(t2_size);
            op_t_[t1_size-1][t2_size-1].resize(s_size);
            op_t_[t1_size-1][t2_size-1][s_size-1].push_back(t);
            t_size = 1;
            op_values_[t1_size-1].resize(t2_size);
            op_values_[t1_size-1][t2_size-1].resize(s_size);
            op_values_[t1_size-1][t2_size-1][s_size-1].resize(t_size);
            op_values_[t1_size-1][t2_size-1][s_size-1][t_size-1].push_back(e);
            target_t2 = t2;
            target_s = s;
            target_t = t;
            continue;
        }

        if (s != target_s) {
            op_s_[t1_size-1][t2_size-1].push_back(s);
            s_size++;
            op_t_[t1_size-1][t2_size-1].resize(s_size);
            op_t_[t1_size-1][t2_size-1][s_size-1].push_back(t);
            t_size = 1;
            op_values_[t1_size-1][t2_size-1].resize(s_size);
            op_values_[t1_size-1][t2_size-1][s_size-1].resize(t_size);
            op_values_[t1_size-1][t2_size-1][s_size-1][t_size-1].push_back(e);
            target_s = s;
            target_t = t;
            continue;
        }

        if (t != target_t) {
            op_t_[t1_size-1][t2_size-1][s_size-1].push_back(t);
            t_size++;
            op_values_[t1_size-1][t2_size-1][s_size-1].resize(t_size);
            op_values_[t1_size-1][t2_size-1][s_size-1][t_size-1].push_back(e);
            target_t = t;
            continue;
        }

        op_values_[t1_size-1][t2_size-1][s_size-1][t_size-1].push_back(e);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "op index construction time: " << (double)duration.count()/1000 << " ms" << std::endl;
    print_peak_memory_usage();
}

void MSTIndex::create_qt_idx() {
    auto start = std::chrono::high_resolution_clock::now();
    //create_bl_idx_plus();
    // (s, t1, t2, t)
    Rectangle boundary(0, -1, time_window_sz, time_window_sz);
    qt_idx_ = std::make_unique<QuadTree>(boundary);
    for (auto i = 0; i < time_window_sz; i++) {
        for (auto j = 0; j < bl_keys_[i].size(); j++) {
            for (auto k = 0; k < bl_values_[i][j].size(); k++) {
                Quadruple item;
                if (k == bl_values_[i][j].size() - 1) {
                    item = Quadruple(bl_values_[i][j][k].first, bl_values_[i][j][k].second, INT_MAX, bl_keys_[i][j]);
                } else {
                    item = Quadruple(bl_values_[i][j][k].first, bl_values_[i][j][k].second, bl_values_[i][j][k+1].second, bl_keys_[i][j]);
                }
    
                if (std::get<3>(item).timestamp != std::get<0>(item)) {
                    // t1, s, t2-t1, t-s
                    Rectangle rect(std::get<1>(item), std::get<0>(item),
                                   std::get<2>(item) - std::get<1>(item),
                                   std::get<3>(item).timestamp - std::get<0>(item));
                    qt_idx_->insert(rect, std::get<3>(item));
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "qt index construction time: " << (double)duration.count()/1000 << " ms" << std::endl;
    print_peak_memory_usage();

    //std::cout << "Quad-tree status: " << qt_idx_->check() << std::endl;
}

void MSTIndex::qt_idx_count(std::string &index_file) {
    create_bl_idx_plus();
    // (s, t1, t2, t)
    std::vector<Quadruple> items;
    for (auto i = 0; i < time_window_sz; i++) {
        for (auto j = 0; j < bl_keys_[i].size(); j++) {
            for (auto k = 0; k < bl_values_[i][j].size(); k++) {
                if (k == bl_values_[i][j].size() - 1) {
                    items.emplace_back(bl_values_[i][j][k].first, bl_values_[i][j][k].second, INT_MAX, bl_keys_[i][j]);
                } else {
                    items.emplace_back(bl_values_[i][j][k].first, bl_values_[i][j][k].second, bl_values_[i][j][k+1].second, bl_keys_[i][j]);
                }
            }
        }
    }

    Rectangle boundary(0, -1, time_window_sz, time_window_sz);
    qt_idx_ = std::make_unique<QuadTree>(boundary);

    std::unordered_map<int, int> freq;
    for (auto &item : items) {
        // s == t
        if (std::get<3>(item).timestamp == std::get<0>(item)) continue;
        // t1, s, t2-t1, t-s
        Rectangle rect(std::get<1>(item), std::get<0>(item),
                       std::get<2>(item) - std::get<1>(item),
                       std::get<3>(item).timestamp - std::get<0>(item));

        auto cnt = qt_idx_->count(rect);
        auto it = freq.find(cnt);
        if (it == freq.end()) {
            freq.insert(std::make_pair(cnt, 1));
        } else {
            it->second++;
        }
    }

    std::ofstream outFile(index_file);

    if (!outFile) {
        std::cerr << "can not open the file" << std::endl;
    }

    for (const auto& p : freq) {
        outFile << p.first << " " << p.second << std::endl;
    }

    outFile.close();
}

void MSTIndex::create_hy_idx(bool xFirst) {
    auto start = std::chrono::high_resolution_clock::now();
    //create_bl_idx_plus();

    Line boundary(0, 0);
    if (xFirst) {
        boundary = Line(0, time_window_sz);
    } else {
        boundary = Line(-1, time_window_sz-1);
    }
    hy_idx_ = std::make_unique<SegmentTree>(boundary, xFirst);

    // (s, t1, t2, t)
    for (auto i = 0; i < time_window_sz; i++) {
        for (auto j = 0; j < bl_keys_[i].size(); j++) {
            for (auto k = 0; k < bl_values_[i][j].size(); k++) {
                Quadruple item;
                if (k == bl_values_[i][j].size() - 1) {
                    item = Quadruple(bl_values_[i][j][k].first, bl_values_[i][j][k].second, INT_MAX, bl_keys_[i][j]);
                } else {
                    item = Quadruple(bl_values_[i][j][k].first, bl_values_[i][j][k].second, bl_values_[i][j][k+1].second, bl_keys_[i][j]);
                }

                if (std::get<3>(item).timestamp != std::get<0>(item)) {
                    // t1, s, t2-t1, t-s
                    Line x(std::get<1>(item), std::get<2>(item));
                    Line y(std::get<0>(item), std::get<3>(item).timestamp);
                    Edge e = std::get<3>(item);
                    hy_idx_->insert(std::move(x), std::move(y), std::move(e));
                }
            }
        }
    }

    hy_idx_->construct();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "hybrid " << xFirst << " index construction time: " << (double)duration.count()/1000 << " ms" << std::endl;
    print_peak_memory_usage();
}

std::vector<std::vector<Edge>> MSTIndex::query_online(std::vector<TimeWindow> &time_windows) const {
    std::cout << "online query..." << std::endl;

    std::vector<std::vector<Edge>> msts;

    int cnt = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (auto window : time_windows) {
        if (cnt % 100 == 0) {
            // std::cout << "the " << cnt << "th query" << std::endl;
        }
        std::vector<Edge> snapshot;
        for (bint i = window.first; i <= window.second; i++) {
            std::vector<Edge> edges = graph_->get_temporal_edges(i);
            snapshot.insert(snapshot.end(), edges.begin(), edges.end());
        }
        Kruskal_MST(snapshot, graph_->get_num_nodes());
        //msts.emplace_back(Kruskal_MST(snapshot, graph_->get_num_nodes()));
        cnt++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Average time: " << (double)duration.count()/1000/static_cast<double>(time_windows.size()) << " ms" << std::endl;

    return msts;
}

std::vector<std::vector<Edge>> MSTIndex::query_online_global(std::vector<TimeWindow> &time_windows) const
{
    std::cout << "online query (global sorted list)..." << std::endl;

    std::vector<std::vector<Edge>> msts;
    int cnt = 0;
    std::vector<Edge> global_list;
    for (int i = 0; i < time_window_sz; i++)
    {
        std::vector<Edge> edges = graph_->get_temporal_edges(i);
        global_list.insert(global_list.end(), edges.begin(), edges.end());
    }

    std::sort(global_list.begin(), global_list.end(), [](const Edge& a, const Edge& b) {
        if (a.weight != b.weight) {
            return a.weight < b.weight;
        } else if (a.timestamp != b.timestamp) {
            return a.timestamp > b.timestamp;
        } else if (a.src_id != b.src_id) {
            return a.src_id < b.src_id;
        } else {
            return a.dest_id < b.dest_id;
        }
    });

    auto start = std::chrono::high_resolution_clock::now();
    for (auto window : time_windows) {
        if (cnt % 100 == 0) {
            // std::cout << "the " << cnt << "th query" << std::endl;
        }
        std::vector<Edge> snapshot;
        for (auto edge : global_list)
        {
            if (edge.timestamp >= window.first && edge.timestamp <= window.second)
            {
                snapshot.push_back(edge);
            }
        }
        Kruskal_MST_wo_sort(snapshot, graph_->get_num_nodes());
        //msts.emplace_back(Kruskal_MST(snapshot, graph_->get_num_nodes()));
        cnt++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Average time: " << (double)duration.count()/1000/static_cast<double>(time_windows.size()) << " ms" << std::endl;

    return msts;
}

std::vector<Edge> MSTIndex::k_way_merge(std::vector<std::vector<Edge>>& sorted_lists)
{
    using PQElement = std::pair<Edge, std::pair<int, int>>; // (edge, (list_idx, edge_idx))
    auto cmp = [](const PQElement& a, const PQElement& b) {
        return a.first.weight > b.first.weight;  // min-heap based on weight
    };

    std::priority_queue<PQElement, std::vector<PQElement>, decltype(cmp)> min_heap(cmp);
    std::vector<Edge> merged_result;

    // Initialize: push the first element of each list
    for (int i = 0; i < sorted_lists.size(); ++i) {
        if (!sorted_lists[i].empty()) {
            min_heap.emplace(sorted_lists[i][0], std::make_pair(i, 0));
        }
    }

    while (!min_heap.empty()) {
        auto [edge, pos] = min_heap.top();
        min_heap.pop();
        merged_result.push_back(edge);

        int list_idx = pos.first;
        int edge_idx = pos.second;

        // Push the next edge from the same list, if any
        if (edge_idx + 1 < sorted_lists[list_idx].size()) {
            min_heap.emplace(sorted_lists[list_idx][edge_idx + 1], std::make_pair(list_idx, edge_idx + 1));
        }
    }

    return merged_result;
}

std::vector<std::vector<Edge>> MSTIndex::query_online_k_way(std::vector<TimeWindow> &time_windows)
{
    std::cout << "online query (k way)..." << std::endl;

    std::vector<std::vector<Edge>> msts;
    std::vector<std::vector<Edge>> basic_lists;
    int cnt = 0;
    for (int i = 0; i < time_window_sz; i++)
    {
        std::vector<Edge> edges = graph_->get_temporal_edges(i);
        std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
            if (a.weight != b.weight) {
                return a.weight < b.weight;
            } else if (a.timestamp != b.timestamp) {
                return a.timestamp > b.timestamp;
            } else if (a.src_id != b.src_id) {
                return a.src_id < b.src_id;
            } else {
                return a.dest_id < b.dest_id;
            }
        });
        basic_lists.push_back(edges);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (auto window : time_windows) {
        if (cnt % 100 == 0) {
            // std::cout << "the " << cnt << "th query" << std::endl;
        }
        std::vector<Edge> snapshot;
        std::vector<std::vector<Edge>> k_lists;
        for (auto i = window.first; i <= window.second; i++)
        {
            k_lists.push_back(basic_lists[i]);
        }
        snapshot = k_way_merge(k_lists);
        Kruskal_MST_wo_sort(snapshot, graph_->get_num_nodes());
        //msts.emplace_back(Kruskal_MST(snapshot, graph_->get_num_nodes()));
        cnt++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Average time: " << (double)duration.count()/1000/static_cast<double>(time_windows.size()) << " ms" << std::endl;

    return msts;
}

std::vector<std::vector<Edge>> MSTIndex::query_bl(std::vector<TimeWindow> &time_windows) {
    auto valid = [](ExpireTimeRecord record, bint t_begin, bint t_end) -> bool {
        if (record.back().second <= t_end) {
            if (record.back().first < t_begin) {
                return true;
            } else {
                return false;
            }
        }
        for (int i = 0; i < record.size() - 1; i++) {
            if (record[i+1].second > t_end) {
                if (record[i].first < t_begin) {
                    return true;
                } else {
                    return false;
                }
            }
        }
        return false;
    };

    std::cout << "query baseline..." << std::endl;
    std::vector<std::vector<Edge>> msts(time_windows.size());
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < time_windows.size(); i++) {
        if (i % 100 == 0) {
            // std::cout << "the " << i << "th query" << std::endl;
        }
        TimeWindow time_window = time_windows[i];
        std::vector<Edge> mst;
        for (auto j = time_window.first; j <= time_window.second; j++) {
            for (auto k = 0; k < bl_keys_[j].size(); k++) {
                Edge edge = bl_keys_[j][k];
                auto record = bl_values_[j][k];
                if (valid(record, time_window.first, time_window.second)) {
                    mst.emplace_back(edge);
                }
            }
        }
        //msts[i] = mst;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Average time: " << (double)duration.count()/1000/static_cast<double>(time_windows.size()) << " ms" << std::endl;

    return msts;
}

std::vector<std::vector<Edge>> MSTIndex::query_op(std::vector<TimeWindow> &time_windows) {
    std::cout << "query optimized..." << std::endl;
    std::vector<std::vector<Edge>> msts(time_windows.size());
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < time_windows.size(); i++) {
        if (i % 100 == 0) {
            //std::cout << "the " << i << "th query" << std::endl;
        }
        TimeWindow time_window = time_windows[i];
        std::vector<Edge> mst;
        auto t_s = time_window.first;
        auto t_e = time_window.second;
        for (auto j = 0; j < op_t1_.size(); j++) {
            auto t1 = op_t1_[j];
            if (t1 > t_e) break;
            for (auto k = 0; k < op_t2_[j].size(); k++) {
                auto t2 = op_t2_[j][k];
                if (t2 <= t_e) break;
                for (auto l = 0; l < op_s_[j][k].size(); l++) {
                    auto s = op_s_[j][k][l];
                    if (t_s <= s) break;
                    for (auto m = 0; m < op_t_[j][k][l].size(); m++) {
                        auto t = op_t_[j][k][l][m];
                        if (t < t_s) break;
                        mst.insert(mst.end(), op_values_[j][k][l][m].begin(), op_values_[j][k][l][m].end());
                    }
                }
            }
        }
        //msts[i] = mst;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Average time: " << (double)duration.count()/1000/static_cast<double>(time_windows.size()) << " ms" << std::endl;

    return msts;
}

std::vector<std::vector<Edge>> MSTIndex::query_qt(std::vector<TimeWindow> &time_windows) {
    std::cout << "query qt..." << std::endl;
    std::vector<std::vector<Edge>> msts(time_windows.size());
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < time_windows.size(); i++) {
        if (i % 100 == 0) {
            // std::cout << "the " << i << "th query" << std::endl;
        }
        TimeWindow time_window = time_windows[i];
        std::vector<Edge> mst;
        auto t_s = time_window.first;
        auto t_e = time_window.second;
        qt_idx_->retrieve(Point(t_e, t_s), mst);
        //msts[i] = mst;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Average time: " << (double)duration.count()/1000/static_cast<double>(time_windows.size()) << " ms" << std::endl;

    return msts;
}

std::vector<std::vector<Edge>> MSTIndex::query_hy(std::vector<TimeWindow> &time_windows) {
    std::cout << "query hybrid..." << std::endl;
    std::vector<std::vector<Edge>> msts(time_windows.size());
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < time_windows.size(); i++) {
        if (i % 100 == 0) {
            // std::cout << "the " << i << "th query" << std::endl;
        }
        TimeWindow time_window = time_windows[i];
        std::vector<Edge> mst;
        auto t_s = time_window.first;
        auto t_e = time_window.second;
        if (hy_idx_->xFirst) {
            hy_idx_->retrieve(Point(t_e, t_s), mst);
        } else {
            hy_idx_->retrieve(Point(t_s, t_e), mst);
        }
        //msts[i] = mst;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Average time: " << (double)duration.count()/1000/static_cast<double>(time_windows.size()) << " ms" << std::endl;

    return msts;
}

void MSTIndex::validate(std::vector<TimeWindow> &time_windows, std::string index, std::string ratio) {
    std::cout << "validate index..." << std::endl;

    std::vector<std::vector<Edge>> index_result;
    std::vector<std::vector<Edge>> ground_truth = query_online(time_windows);
    int cnt = 0;

    while (cnt < 4) {
        if (cnt == 0) {
            std::string bl_index = index+"_bl_"+ratio+".bin";
            load_bl_idx(bl_index);
            index_result = query_bl(time_windows);
        } else if (cnt == 1) {
            std::string op_index = index+"_op_"+ratio+".bin";
            load_op_idx(op_index);
            index_result = query_op(time_windows);
        } else if (cnt == 2) {
            std::string qt_index = index+"_qt_"+ratio+".bin";
            load_qt_idx(qt_index);
            index_result = query_qt(time_windows);
        } else if (cnt == 3) {
            std::string hy_index = index+"_hy1_"+ratio+".bin";
            load_hy_idx(hy_index);
            index_result = query_hy(time_windows);
        }

        for (auto i = 0; i < time_windows.size(); i++) {
            std::set<Edge> ground_truth_set;
            for (auto &e : ground_truth[i]) {
                ground_truth_set.insert(e);
            }
            std::set<Edge> index_result_set;
            for (auto &e : index_result[i]) {
                index_result_set.insert(e);
            }
            if (ground_truth_set == index_result_set) {
                continue;
            } else {
                int tf = 0;
                int ft = 0;
                for (auto& e : ground_truth_set) {
                    if (index_result_set.find(e) == index_result_set.end()) {
                        tf++;
                        std::cout << "true answer " << e << " not in [" << time_windows[i].first << ", " << time_windows[i].second << "]" << std::endl;
                    }
                }
                for (auto& e : index_result_set) {
                    if (ground_truth_set.find(e) == ground_truth_set.end()) {
                        ft++;
                        std::cout << "false answer " << e << " in [" << time_windows[i].first << ", " << time_windows[i].second << "]" << std::endl;
                    }
                }
                // std::cout << tf << " " << ft << std::endl;
                if (tf + ft > 0) {
                    exit(-1);
                }
            }
        }
        if (cnt == 1) {
            delete_op_idx();
        } else if (cnt == 2) {
            delete_qt_idx();
        } else if (cnt == 3) {
            delete_hy_idx();
        }
        cnt++;
    }
}

void MSTIndex::save_bl_idx(std::string &index_file) {
    std::cout << "save bl index..." << std::endl;
    std::ofstream file(index_file, std::ios::binary);

    if (!file) {
        throw std::runtime_error("can not open given file");
    }

    // bl_keys_
    save_nested_vector(file, bl_keys_);
    // bl_values_
    save_double_nested_vector(file, bl_values_);

    file.close();
    std::ifstream size_reader(index_file, std::ios::binary | std::ios::ate);
    if (!size_reader) {
        std::cout << "Failed to open file: " << index_file << std::endl;
        exit(-1);
    }
    std::streampos index_size = size_reader.tellg();
    auto m_size = index_size/1024/1024;
    std::cout << "bl index size: " << m_size << "MB" << std::endl;
}

void MSTIndex::load_bl_idx(std::string &index_file) {
    std::cout << "load bl index..." << std::endl;
    std::ifstream file(index_file, std::ios::binary);

    if (!file) {
        throw std::runtime_error("can not open index file");
    }

    // bl_keys_
    load_nested_vector(file, bl_keys_);
    // bl_values_
    load_double_nested_vector(file, bl_values_);
    time_window_sz = 2000;

    auto s_size = 0;
    for (auto i = 0; i < bl_keys_.size(); i++)
    {
        s_size += bl_keys_[i].size();
    }

    size_t t_size = 0;
    for (auto i = 0; i < bl_values_.size(); i++) {
        for (auto j = 0; j < bl_values_[i].size(); j++)
        {
            auto bl_value = bl_values_[i][j];
            if (bl_value.back().first != i)
            {
                t_size += time_window_sz - i;
            } else
            {
                t_size += bl_value.back().second - i;
            }
        }
    }

    std::cout << "t star is: " <<  (double)t_size/s_size << std::endl;

    file.close();
}

template <typename T>
void MSTIndex::save_vector(std::ofstream &out, const std::vector<T> &vec) {
    size_t size = vec.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(T));
}

template <typename T>
void MSTIndex::save_nested_vector(std::ofstream &out, const std::vector<std::vector<T>> &nested_vec) {
    size_t size = nested_vec.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const std::vector<T> &vec : nested_vec) {
        save_vector(out, vec);
    }
}

template <typename T>
void MSTIndex::save_double_nested_vector(std::ofstream &out, const std::vector<std::vector<std::vector<T>>> &nested_vec) {
    size_t size = nested_vec.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const std::vector<std::vector<T>> &vec : nested_vec) {
        save_nested_vector(out, vec);
    }
}

template <typename T>
void MSTIndex::save_triple_nested_vector(std::ofstream &out, const std::vector<std::vector<std::vector<std::vector<T>>>> &nested_vec) {
    size_t size = nested_vec.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const std::vector<std::vector<std::vector<T>>> &vec : nested_vec) {
        save_double_nested_vector(out, vec);
    }
}

template <typename T>
void MSTIndex::save_quadruple_nested_vector(std::ofstream &out, const std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>> &nested_vec) {
    size_t size = nested_vec.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const std::vector<std::vector<std::vector<std::vector<T>>>> &vec : nested_vec) {
        save_triple_nested_vector(out, vec);
    }
}

void MSTIndex::save_op_idx(std::string &index_file) {
    std::cout << "save op index..." << std::endl;
    std::ofstream file(index_file, std::ios::binary);

    if (!file) {
        throw std::runtime_error("can not open given file");
    }

    // op_t1_
    save_vector(file, op_t1_);
    // op_t2_
    save_nested_vector(file, op_t2_);
    // op_s_
    save_double_nested_vector(file, op_s_);
    // op_t_
    save_triple_nested_vector(file, op_t_);
    // op_values_
    save_quadruple_nested_vector(file, op_values_);

    file.close();
    std::ifstream size_reader(index_file, std::ios::binary | std::ios::ate);
    if (!size_reader) {
        std::cout << "Failed to open file: " << index_file << std::endl;
        exit(-1);
    }
    std::streampos index_size = size_reader.tellg();
    auto m_size = index_size/1024/1024;
    std::cout << "op index size: " << m_size << "MB" << std::endl;

    // release memory
    delete_op_idx();
}

void MSTIndex::save_qt_idx(std::string &index_file) {
    std::cout << "save qt index..." << std::endl;

    qt_idx_->serialize(index_file);

    std::ifstream size_reader(index_file, std::ios::binary | std::ios::ate);
    if (!size_reader) {
        std::cout << "Failed to open file: " << index_file << std::endl;
        exit(-1);
    }
    std::streampos index_size = size_reader.tellg();
    auto m_size = index_size/1024/1024;
    std::cout << "qt index size: " << m_size << "MB" << std::endl;

    // release memory
    delete_qt_idx();
}

void MSTIndex::save_hy_idx(std::string &index_file) {
    std::cout << "save hybrid index..." << std::endl;

    hy_idx_->serialize(index_file);

    std::ifstream size_reader(index_file, std::ios::binary | std::ios::ate);
    if (!size_reader) {
        std::cout << "Failed to open file: " << index_file << std::endl;
        exit(-1);
    }
    std::streampos index_size = size_reader.tellg();
    auto m_size = index_size/1024/1024;
    std::cout << "hybrid " << hy_idx_->xFirst << " index size: " << m_size << "MB" << std::endl;

    // release memory
    delete_hy_idx();
}

template <typename T>
void MSTIndex::load_vector(std::ifstream &in, std::vector<T> &vec) {
    size_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);
    in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(T));
}

template <typename T>
void MSTIndex::load_nested_vector(std::ifstream &in, std::vector<std::vector<T>> &nested_vec) {
    size_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    nested_vec.resize(size);
    for (std::vector<T> &vec : nested_vec) {
        load_vector(in, vec);
    }
}

template <typename T>
void MSTIndex::load_double_nested_vector(std::ifstream &in, std::vector<std::vector<std::vector<T>>> &nested_vec) {
    size_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    nested_vec.resize(size);
    for (std::vector<std::vector<T>> &vec : nested_vec) {
        load_nested_vector(in, vec);
    }
}

template <typename T>
void MSTIndex::load_triple_nested_vector(std::ifstream &in, std::vector<std::vector<std::vector<std::vector<T>>>> &nested_vec) {
    size_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    nested_vec.resize(size);
    for (std::vector<std::vector<std::vector<T>>> &vec : nested_vec) {
        load_double_nested_vector(in, vec);
    }
}

template <typename T>
void MSTIndex::load_quadruple_nested_vector(std::ifstream &in, std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>> &nested_vec) {
    size_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    nested_vec.resize(size);
    for (std::vector<std::vector<std::vector<std::vector<T>>>> &vec : nested_vec) {
        load_triple_nested_vector(in, vec);
    }
}

void MSTIndex::load_op_idx(std::string &index_file) {
    std::cout << "load op index..." << std::endl;
    std::ifstream file(index_file, std::ios::binary);

    if (!file) {
        throw std::runtime_error("can not open index file");
    }

    // op_t1_
    load_vector(file, op_t1_);
    // op_t2_
    load_nested_vector(file, op_t2_);
    // op_s_
    load_double_nested_vector(file, op_s_);
    // op_t_
    load_triple_nested_vector(file, op_t_);
    // op_values_
    load_quadruple_nested_vector(file, op_values_);

    auto s_size = 0;
    for (auto i = 0; i < op_s_.size(); i++) {
        for (auto j = 0; j < op_s_[i].size(); j++) {
            s_size += op_s_[i][j].size();
        }
    }

    auto t_size = 0;
    for (auto i = 0; i < op_t_.size(); i++) {
        for (auto j = 0; j < op_t_[i].size(); j++) {
            for (auto k = 0; k < op_t_[i][j].size(); k++) {
                t_size += op_t_[i][j][k].size();
            }
        }
    }

    std::cout << "t bar is: " <<  (double)t_size/s_size << std::endl;

    file.close();
}

void MSTIndex::load_qt_idx(std::string &index_file) {
    std::cout << "load qt index..." << std::endl;

    qt_idx_ = std::make_unique<QuadTree>(Rectangle());
    qt_idx_->deserialize(index_file);
}

void  MSTIndex::load_hy_idx(std::string &index_file) {
    std::cout << "load hy index..." << std::endl;
    hy_idx_ = std::make_unique<SegmentTree>(Line(0, 0), true);
    hy_idx_->deserialize(index_file);
}

void MSTIndex::delete_op_idx() {
    op_t1_.clear();
    op_t1_.shrink_to_fit();

    op_t2_.clear();
    op_t2_.shrink_to_fit();

    op_s_.clear();
    op_s_.shrink_to_fit();

    op_t_.clear();
    op_t_.shrink_to_fit();

    op_values_.clear();
    op_values_.shrink_to_fit();
}

void MSTIndex::delete_qt_idx() {
    qt_idx_->clear();
}

void MSTIndex::delete_hy_idx() {
    hy_idx_->clear();
}