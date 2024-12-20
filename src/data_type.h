#ifndef DATA_TYPE_H
#define DATA_TYPE_H

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <cstring>

typedef int bint;
typedef std::pair<bint, bint> TimeWindow;
typedef std::vector<std::pair<bint, bint>> ExpireTimeRecord;

/*
 * Edge
 */
struct Edge {
    int src_id = 0;
    int dest_id = 0;
    int weight = 0;
    bint timestamp = 0;
    Edge(int src_id, int dest_id, int weight, bint timestamp) : src_id(src_id), dest_id(dest_id), weight(weight), timestamp(timestamp) {}
    Edge() = default;
    // hashmap
    bool operator==(const Edge& other) const {
        return (src_id == other.src_id && dest_id == other.dest_id && weight == other.weight && timestamp == other.timestamp);
    }
    // set
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Edge& a) {
    os << "(" << a.src_id << ", " << a.dest_id << ", " << a.weight << ", " << a.timestamp << ")";
    return os;
}

struct edge_hash {
    std::size_t operator()(const Edge& edge) const {
        return std::hash<float>()(edge.weight);
    }
};

typedef std::tuple<bint, bint, bint, Edge> Quadruple;

inline bool compareQuadruples(const Quadruple& a, const Quadruple& b) {
    // t1 from small to large
    if (std::get<1>(a) != std::get<1>(b)) {
        return std::get<1>(a) < std::get<1>(b);
    }
    // t2 from large to small
    if (std::get<2>(a) != std::get<2>(b)) {
        return std::get<2>(a) > std::get<2>(b);
    }
    // s from small to large
    if (std::get<0>(a) != std::get<0>(b)) {
        return std::get<0>(a) < std::get<0>(b);
    }
    // t from large to small
    return std::get<3>(a).timestamp > std::get<3>(b).timestamp;
}

class UnionFind {
public:
    explicit UnionFind(int n) {
        parent.reserve(n);
        parent.resize(n);
        rank.reserve(n);
    }

    void clear(int n) {
        rank.assign(n, 0);
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    int find(int x) {
        int root = x;
        // find root node
        while (root != parent[root]) {
            root = parent[root];
        }
        // path compression
        while (x != root) {
            int next = parent[x];
            parent[x] = root;
            x = next;
        }
        return root;
    }

    void unite(int rootX, int rootY) {
        if (rootX != rootY) {
            if (rank[rootX] > rank[rootY]) {
                parent[rootY] = rootX;
            } else if (rank[rootX] < rank[rootY]) {
                parent[rootX] = rootY;
            } else {
                parent[rootY] = rootX;
                rank[rootX] += 1;
            }
        }
    }

    bool connected(int x, int y, int &rootX, int &rootY) {
        rootX = find(x);
        rootY = find(y);
        return rootX == rootY;
    }

private:
    std::vector<int> parent;
    std::vector<int> rank;
};

#endif // DATA_TYPE_H