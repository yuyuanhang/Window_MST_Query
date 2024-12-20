#ifndef DYNAMICMST_H_
#define DYNAMICMST_H_

#include <string>
#include <vector>
#include <cstring>
#include <climits>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <cmath>
#include <set>
#include <queue>
#include <iostream>

#define MAXST 1024
#define INF 2000000000

using namespace std;

namespace DMST {
    inline long candsize = 0;
    inline long candtimes = 0;
    inline long nontreesize = 0;
    inline long nontreetimes = 0;
    inline long qqsize = 0;
    inline long qqtimes = 0;
    inline long einf = 0;
    inline long qcheck = 0;
    inline long aa = 0;
    inline long bb = 0;
    inline long cc = 0;
    inline long dd = 0;

    class Edge {
    public:
        int u, v, w;

    public:
        Edge();

        Edge(int u, int v, int w);
    };

    struct EdgeHash {
        std::size_t operator()(const DMST::Edge& e) const {
            std::size_t h1 = std::hash<int>{}(std::min(e.u, e.v));
            std::size_t h2 = std::hash<int>{}(std::max(e.u, e.v));
            std::size_t h3 = std::hash<int>{}(e.w);
            return h1 ^ (h2 << 1) ^ (h3 << 2);  // Combine hashes
        }
    };

    typedef pair<Edge, int> Pe;

    class Adj {
    public:
        int v, w;
    public:
        Adj();

        Adj(int v, int w);
    };

    class Node {
    public:
        //for graph
        vector<Adj> adj;

        //for tree
        Edge p; //parent node in the tree, -1 means the node is a root node
        int sub_cnt; //number of descendants in the tree

        //replacement edge
        Edge min_e;

    public:
        Node();

    public:
        vector<Adj> del_buf;
        vector<Adj> ins_buf;

    public:
        void insert_adj(int v, int w);

        void delete_adj(int v, int w);

        void flush();
    };

    class DynamicMST {
    public:
        static void create_bin(string path, bool multi_edge = false);

        static bool get_edge(char *line, int &a, int &b, int &w, int num_cnt);

        static bool get_temporal_edge(char *line, int &a, int &b, int &t, int num_cnt);

        static int get_num_cnt(string path);

        static void create_stream(string path, bool is_bipartite);

    public:
        int n, rst;
        vector<int> cand, l, status, cand1, cand2, c;
        vector<Node> nodes;

        vector<bool> used;

    public:
        long long s1, s2, s3, t, t3;

    public:
        DynamicMST(string path, bool load_graph = true);

        DynamicMST(int node_count);

        void init();

        int update_min_e(int u, int v, int w);

        void search_min_e(vector<int> &cand, int u, Edge e_old);

        int check_status(int u);

        void search_min_e(vector<int> &cand, Edge e_old);


        void reroot(int u);

        void _reroot(int u);

        bool insert_edge_in_graph(int u, int v, int w);

        bool delete_edge_in_graph(int u, int v, int w);

        int insert_edge(int u, int v, int w, Edge &removed_e);

        int delete_edge(int u, int v, int w);

        int _insert_edge(int u, int v, int w, Edge &removed_e);

        int _delete_edge(int u, int v, int w);

        void replace(int u, int v, int w, int _u, int _v, int _w);

        void sample_edges(vector<Edge> &edges, int cnt);

        long long mst_weight();

        void print_mst();

        bool is_decedent(int u, int v);

        bool is_out_edge(int u, Edge e);
    };

//=================== implementation ==================

    inline int get_f(vector<int> &f, int v) {
        if (f[v] != v) f[v] = get_f(f, f[v]);
        return f[v];
    }

    inline void union_f(vector<int> &f, int a, int b) {
        int fa = get_f(f, a);
        int fb = get_f(f, b);
        f[fa] = fb;
    }

    inline bool operator<(const Edge &x, const Edge &y) {
        if (x.w < y.w) return true;
        if (x.w > y.w) return false;
        int minx = min(x.u, x.v), miny = min(y.u, y.v);
        if (minx < miny) return true;
        if (minx > miny) return false;
        return max(x.u, x.v) < max(y.u, y.v);
    }


    inline bool operator==(const Edge &x, const Edge &y) {
        if (x.w != y.w) return false;
        int minx = min(x.u, x.v), miny = min(y.u, y.v);
        if (minx != miny) return false;
        return max(x.u, x.v) == max(y.u, y.v);
    }

    inline bool operator<(const Adj &x, const Adj &y) {
        if (x.w < y.w) return true;
        if (x.w > y.w) return false;
        return x.v < y.v;
    }


    inline bool operator==(const Adj &x, const Adj &y) {
        return x.w == y.w && x.v == y.v;
    }
}

#endif /* DYNAMICMST_H_ */