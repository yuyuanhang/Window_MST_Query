//
// Created by yuanyu on 28/08/24.
//
#include "dynamic_mst.h"

namespace DMST {
    Edge::Edge() {
        u = -1;
        v = -1;
        w = INF;
    }

    Edge::Edge(int u, int v, int w) {
        this->u = u;
        this->v = v;
        this->w = w;
    }

    Adj::Adj() {v = -1; w = INF;}
    Adj::Adj(int v, int w) {this->v = v; this->w = w;}


    Node::Node() {
        p = Edge(); sub_cnt = 0;
    }


    void Node::insert_adj(int v, int w) {
        ins_buf.push_back(Adj(v,w));
    }
    void Node::delete_adj(int v, int w) {
        del_buf.push_back(Adj(v,w));
    }
    void Node::flush() {
        if(ins_buf.size() == 0 && del_buf.size() == 0) return;
        sort(ins_buf.begin(), ins_buf.end());
        sort(del_buf.begin(), del_buf.end());

        int i = 0, d = 0, ni = (int)ins_buf.size(), nd = (int)del_buf.size();

        vector<Adj> l;
        for(int j = 0; j <= (int) adj.size(); ++j) {
            Adj v = j < (int) adj.size() ? adj[j] : Adj();
            for(;i < ni && ins_buf[i] < v; ++i) {
                if(d < nd && del_buf[d] == ins_buf[i]) {++d; continue;}
                l.push_back(ins_buf[i]);
            }
            if(j == (int) adj.size()) break;
            if(d < nd && del_buf[d] == v) {++d; continue;}
            l.push_back(v);
        }
        ins_buf.clear(); del_buf.clear();
        adj = l;
    }

    bool DynamicMST::insert_edge_in_graph(int u, int v, int w) {
        if(u < 0 || u >= n || v < 0 || v >= n || u == v) return false;
        nodes[u].insert_adj(v, w);
        nodes[v].insert_adj(u, w);
        return true;
    }

    bool DynamicMST::delete_edge_in_graph(int u, int v, int w) {
        if(u < 0 || u >= n || v < 0 || v >= n || u == v) return false;
        nodes[u].delete_adj(v, w);
        nodes[v].delete_adj(u, w);
        return true;
    }

    int DynamicMST::insert_edge(int u, int v, int w, Edge &removed_e) {
        if(!insert_edge_in_graph(u,v,w)) return -1;
        return _insert_edge(u,v,w, removed_e);
    }

    int DynamicMST::delete_edge(int u, int v, int w) {
        if(!delete_edge_in_graph(u,v,w)) return -1;
        return _delete_edge(u,v,w);
    }

    void DynamicMST::replace(int u, int v, int w, int _u, int _v, int _w) {
        reroot(_u);
        if(nodes[u].p.v != v) swap(u,v);
        for(int z = v; z != -1; z = nodes[z].p.v) nodes[z].sub_cnt -= nodes[u].sub_cnt;
        nodes[u].p = Edge();
        reroot(_v);

        if(nodes[_v].sub_cnt > nodes[_u].sub_cnt){
            swap(_u, _v);
        }
        nodes[_v].p = Edge(_v,_u,_w);
        nodes[_u].sub_cnt += nodes[_v].sub_cnt;

        update_min_e(u,v,w);
    }


    int DynamicMST::_insert_edge(int u, int v, int w, Edge &removed_e) {
        Edge max_e = Edge(u,v,w);
        int fu = u, fv = v;
        while(fu>=0 && fv>=0 && fu!=fv)
            if(nodes[fu].sub_cnt < nodes[fv].sub_cnt) {
                if(max_e < nodes[fu].p) max_e = nodes[fu].p;
                fu = nodes[fu].p.v;
            } else {
                if(max_e < nodes[fv].p) max_e = nodes[fv].p;
                fv = nodes[fv].p.v;
            }

        if(fu == -1 || fv == -1) {
            reroot(u);
            reroot(v);
            if(nodes[u].sub_cnt > nodes[v].sub_cnt) swap(u,v);
            nodes[u].p = Edge(u,v,w);
            nodes[v].sub_cnt += nodes[u].sub_cnt;
            return 0;
        } else if(max_e == Edge(u,v,w)) {
            update_min_e(u,v,w);
            return 1;
        }

        replace(max_e.u,max_e.v,max_e.w,u,v,w);
        removed_e = max_e;
        return 2;
    }

    int DynamicMST::_delete_edge(int u, int v, int w) {
        Edge e = Edge(u,v,w);
        rst = 1;
        int sc = 0;
        if(nodes[u].p == e || nodes[v].p == e) {
            if(nodes[v].p == e) swap(u,v);
            if(nodes[u].min_e.v == -1) {
                for(int x = v; x != -1; x = nodes[x].p.v) nodes[x].sub_cnt -= nodes[u].sub_cnt;
                nodes[u].p = Edge();
                return 0;
            }

            sc = nodes[u].sub_cnt;
            replace(u,v,w,nodes[u].min_e.u,nodes[u].min_e.v,nodes[u].min_e.w);
            rst = 2;
        }

        if(nodes[u].sub_cnt>nodes[v].sub_cnt) swap(u,v);

        reroot(v);
        int t_size = nodes[v].sub_cnt;
        if(rst == 2) {
            sc = min(sc, t_size-sc);
            ++t3; s3 += sc;
        }

        cand.clear(); cand1.clear(); cand2.clear();
        for(int fu = u; fu != v; fu = nodes[fu].p.v)
            if(nodes[fu].min_e == e) {
                nodes[fu].min_e = Edge();
                cand.push_back(fu);
            }

        candsize += (int)cand.size();
        candtimes++;

        for(int i = (int) cand.size()-1; i >= 0; --i) {
            int x = cand[i];
            if(nodes[x].sub_cnt * 2 < t_size) {
                for(int j = 0; j <= i; ++j) cand1.push_back(cand[j]);
                break;
            } else {
                int fx = nodes[x].p.v;
                reroot(x);
                cand2.push_back(fx);
            }
        }

        if(cand1.size()) search_min_e(cand1, e);
        if(cand2.size()) search_min_e(cand2, e);

        return rst;
    }

    int DynamicMST::check_status(int u) {
        int rst = INF;
        for(int x = u; x != -1; x = nodes[x].p.v)
            if(status[x] != -1) {rst = status[x]; break;}
        for(int x = u; x != -1; x = nodes[x].p.v) {
            if(status[x] != -1) break;
            status[x] = rst;
            l.push_back(x);
        }

        return rst;
    }

    void DynamicMST::search_min_e(vector<int> &cand, Edge e_old) {
        l.clear();
        for(int i = 0; i < (int) cand.size(); ++i) {
            status[cand[i]] = i;
            l.push_back(cand[i]);
        }

        for(int i = 0; i < (int) cand.size(); ++i)
            search_min_e(cand, cand[i], e_old);

        for(int i = 0; i < (int) l.size(); ++i) status[l[i]] = -1;

        if(rst == 2) {s1 += l.size(); s2 += nodes[cand[cand.size()-1]].sub_cnt; ++t;}
    }

    void DynamicMST::search_min_e(vector<int> &cand, int u, Edge e_old) {
        priority_queue<Pe, vector<Pe>, greater<Pe> > h;
        h.push(make_pair(nodes[u].min_e,u));

        bool qqflag = true;


        while(!h.empty()) {
            int x = h.top().second;
            Edge e = h.top().first;
            h.pop();
            if(qqflag){
                qqtimes++;
                qqflag = false;
            }
            if(e.w==INF)
                einf++;
            qqsize++;
            if(nodes[u].min_e < e) break;
            if(x != u && nodes[x].min_e.v >= 0 )  {
                dd++;
                if (nodes[x].min_e == e)
                    aa++;
                if(!(e<e_old))
                    bb++;
                if(nodes[x].min_e == e && !(e<e_old))
                    cc++;
            }
            if(x != u && nodes[x].min_e.v >= 0 && !(e<e_old))  {
                qcheck++;
                int p = check_status(e.v);
                if(p > status[u]) {
                    for(p=(p==INF?(cand.size()-1):(p-1)); p >= status[u]; --p)
                        if(e < nodes[cand[p]].min_e) nodes[cand[p]].min_e = e;
                    break;
                }
            }

            nodes[x].flush();
            c.clear();
            bool search_adj = true;
            nontreetimes++;
            for (int j = 0; j < (int)nodes[x].adj.size(); ++j)
            {
                int y = nodes[x].adj[j].v;
                if(y == nodes[x].p.v && !used[y]) {used[y] = true; c.push_back(y); continue;}
                if(nodes[y].p.v == x && !used[y]) {
                    used[y] = true;
                    c.push_back(y);
                    if(nodes[y].min_e < nodes[u].min_e)
                        h.push(make_pair(nodes[y].min_e,y));
                    continue;
                }
                if(!search_adj)
                    continue;
                Edge _e = Edge(x, y, nodes[x].adj[j].w);
                if(_e < nodes[u].min_e && (x == u || e < _e) && !(_e<e_old) ) {
                    nontreesize++;
                    int p = check_status(y);
                    if(p > status[u]){
                        for(p=(p==INF?(cand.size()-1):(p-1)); p >= status[u]; --p)
                            if(_e < nodes[cand[p]].min_e) nodes[cand[p]].min_e = _e;
                        search_adj = false;
                    }
                }
            }
            for(int j = 0; j < (int) c.size(); ++j) used[c[j]] = false;
        }
    }

    int DynamicMST::update_min_e(int u, int v, int w) {
        Edge eu = Edge(u,v,w), ev = Edge(v,u,w);
        int cnt = 0;
        for(int fu = u, fv = v; fu != fv;)
            if(nodes[fu].sub_cnt < nodes[fv].sub_cnt) {
                if(eu < nodes[fu].min_e) {nodes[fu].min_e = eu;++cnt;}
                fu = nodes[fu].p.v;
            } else {
                if(ev < nodes[fv].min_e) {nodes[fv].min_e = ev;++cnt;}
                fv = nodes[fv].p.v;
            }
        return cnt;
    }

    void DynamicMST::_reroot(int u) {
        int p = nodes[u].p.v;
        if(p == -1) return;
        _reroot(p);
        nodes[p].p = Edge(p,u,nodes[u].p.w);
        int c = nodes[u].sub_cnt;
        nodes[u].sub_cnt = nodes[p].sub_cnt;
        nodes[p].sub_cnt -= c;
        nodes[p].min_e = Edge(nodes[u].min_e.v, nodes[u].min_e.u, nodes[u].min_e.w);
    }

    void DynamicMST::reroot(int u) {
        _reroot(u);
        nodes[u].p = Edge();
        nodes[u].min_e = Edge();
    }

    void DynamicMST::init() {
        vector<int> f(n);
        vector<Edge> e_all;
        for(int u = 0; u < n; ++u) {
            nodes[u].min_e = Edge();
            nodes[u].p = Edge();
            nodes[u].sub_cnt = 1;
            nodes[u].flush();
            for(int i = 0; i < (int) nodes[u].adj.size(); ++i) {
                int v = nodes[u].adj[i].v, w = nodes[u].adj[i].w;
                if(v > u) e_all.push_back(Edge(u, v, w));
            }
            f[u] = u;
        }

        vector<vector<Adj> > mst(n);

        sort(e_all.begin(), e_all.end());

        int n_comp = n;
        for(size_t i = 0; i < e_all.size(); ++i) {
            int u = e_all[i].u, v = e_all[i].v, w = e_all[i].w;
            int fu = get_f(f,u), fv = get_f(f,v);
            if(fu == fv) continue;
            union_f(f, u, v);
            --n_comp;
            mst[u].push_back(Adj(v,w));
            mst[v].push_back(Adj(u,w));
            e_all[i].u = -1;
        }
        // printf( "n_comp=%d\n", n_comp );

        vector<pair<int,int> > s;
        vector<bool> used(n, false);
        for(int u=0; u<n; ++u) s.push_back(make_pair((int)mst[u].size(),-u));
        sort(s.begin(),s.end());
        vector<int> q;

        for(int i=n-1; i>=0; --i) {
            int f = -s[i].second;

            if(used[f]) continue;
            q.clear();

            used[f] = true;
            q.push_back(f);

            for(int s=0; s < (int)q.size(); ++s) {
                int p = q[s];
                for(int j = 0; j < (int) mst[p].size(); ++j) {
                    int v = mst[p][j].v;
                    if(!used[v]) {
                        used[v] = true;
                        q.push_back(v);
                        nodes[v].p = Edge(v,p,mst[p][j].w);
                    }
                }
            }

            for(int i = (int)q.size()-1; i>0; --i)
                nodes[nodes[q[i]].p.v].sub_cnt += nodes[q[i]].sub_cnt;
        }

        // printf( "updating min_e...\n" );
        for(size_t i = 0; i < e_all.size(); ++i) {
            if((i+1)%10000000 == 0) printf( "%lu/%lu\n", i+1, e_all.size() );
            int u = e_all[i].u, v = e_all[i].v, w = e_all[i].w;
            if(u == -1) continue; //tree edge
            update_min_e(u,v,w);
        }

        // printf( "MST weight = %lld\n", mst_weight() );
    }

    DynamicMST::DynamicMST(string path, bool load_graph) {
        printf( "path = %s\n", path.c_str());

        FILE* fin = fopen( (path+"dgraph.bin").c_str(), "rb" );
        if(fin == NULL || !load_graph) {
            if(fin) fclose(fin);
            fin = fopen( (path+"graph.stream").c_str(), "rb" );
            fread(&n, sizeof(int), 1, fin);
            nodes.resize(n);
            fclose(fin);
        } else {
            fread( &n, sizeof(int), 1, fin );
            nodes.resize(n);
            int *deg = new int[n];
            Adj *dat = new Adj[n];

            printf( "Loading graph...\n" );
            long long m = 0;
            fread( deg, sizeof(int), n, fin );
            for(int i = 0; i < n; ++i) {
                fread(dat, sizeof(Adj), deg[i], fin);
                m += deg[i];
                nodes[i].adj.assign(dat, dat+deg[i]);
            }

            delete[] deg; delete[] dat;
            fclose(fin);
            printf( "Graph loaded, n = %d, m = %lld\n", n, m/2 );
        }

        status.resize(n,-1);
        used.resize(n,false);

        s1 = 0; s2 = 0; s3 = 0; t = 0; t3 = 0; rst = 0;
    }

    DynamicMST::DynamicMST(int node_count) {
        n = node_count;
        nodes.resize(n);

        status.resize(n, -1);
        used.resize(n, false);

        s1 = 0; s2 = 0; s3 = 0; t = 0; t3 = 0; rst = 0;
    }

    bool DynamicMST::get_edge(char *line, int &a, int &b, int &w, int num_cnt) {
        if( !isdigit(line[0]) ) return false;
        vector<char*> v_num;
        int len = (int) strlen(line);
        for( int i = 0; i < len; ++i )
            if( !isdigit(line[i]) && line[i] != '.') line[i] = '\0';
            else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
        if( (int) v_num.size() != num_cnt ) return false;
        sscanf( v_num[0], "%d", &a );
        sscanf( v_num[1], "%d", &b );
        sscanf( v_num[2], "%d", &w );
        return true;
    }

    int DynamicMST::get_num_cnt(string path) {
        FILE *fin = fopen( (path).c_str(), "r" );
        char line[MAXST];
        int cnt = 0, min_cnt = 100;

        while(fgets( line, MAXST, fin ) && cnt < 10) {
            if(!isdigit(line[0])) continue;
            vector<char*> v_num;
            int len = (int) strlen(line);
            for(int i = 0; i < len; ++i)
                if(!isdigit(line[i]) && line[i] != '.') line[i] = '\0';
                else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
            if( (int) v_num.size() < 2 ) continue;
            min_cnt = min(min_cnt, (int) v_num.size());
            ++cnt;
        }
        fclose(fin);
        return min_cnt;
    }

    void DynamicMST::create_bin(string path, bool multi_edge) {
        FILE *fin = fopen((path + "graph.txt").c_str(), "r" );
        char line[MAXST];
        int n = 0, a, b, w, num_cnt = get_num_cnt(path + "graph.txt");
        vector<Edge> el;
        long long cnt = 0, m = 0;

        printf("Loading text, num_cnt = %d, multi_edge = %d ...\n", num_cnt, multi_edge);
        while(fgets( line, MAXST, fin)) {
            if(!get_edge(line, a, b, w, num_cnt)) continue;
            if(a < 0 || b < 0 || a == b) continue;
            el.push_back(Edge(a, b, w));
            n = max(max(n, a+1), b+1);
            if((++cnt) % (long long) 10000000 == 0) printf("%lld lines finished\n", cnt);
        }
        fclose(fin);

        vector<Adj> *con = new vector<Adj>[n];

        for(size_t i = 0; i < el.size(); ++i) {
            con[el[i].u].push_back(Adj(el[i].v, el[i].w));
            con[el[i].v].push_back(Adj(el[i].u, el[i].w));
        }

        vector<bool> used(n,false);
        for(int i = 0; i < n; ++i) {
            sort(con[i].begin(), con[i].end());

            if(!multi_edge) {
                int l = 0;
                for(int j = 0; j < (int) con[i].size(); ++j)
                    if(!used[con[i][j].v]) { con[i][l++] = con[i][j]; used[con[i][j].v] = true;}
                con[i].resize(l);
                for(int j = 0; j < (int) con[i].size(); ++j) used[con[i][j].v] = false;
            }

            m += con[i].size();
        }

        int *deg = new int[n];
        for(int i = 0; i < n; ++i) deg[i] = (int) con[i].size();

        printf( "Saving binary...\n" );

        FILE *fout = fopen((path + "dgraph.bin").c_str(), "wb");
        fwrite(&n, sizeof(int), 1, fout);
        fwrite(deg, sizeof(int), n, fout);
        for(int i = 0; i < n; ++i)
            fwrite(con[i].data(), sizeof(Adj), deg[i], fout);

        fclose(fout);
        printf("Created binary file, n = %d, m = %lld\n", n, m/2);

        delete[] deg; delete[] con;
    }


    void DynamicMST::sample_edges(vector<Edge> &edges, int cnt) {
        long long m = 0;
        for(int v = 0; v < n; ++v) m += nodes[v].adj.size();
        long long now_cnt = 0;
        long long dlt = max((m/2)/cnt,(long long)1);
        edges.clear();
        for(int u = 0; u < n && (int) edges.size() < cnt; ++u) {
            for(int i = 0; i < (int) nodes[u].adj.size(); ++i) {
                int v = nodes[u].adj[i].v;
                if(v < u) continue;
                if(++now_cnt == dlt && (int) edges.size() < cnt) {
                    now_cnt = 0;
                    edges.push_back(Edge(u,v,nodes[u].adj[i].w));
                    if(edges.size() >= cnt) break;
                }
            }
        }
        int s = (int) edges.size();
        for(int i = 0; i < s; ++i) swap(edges[i], edges[rand()%s]);
    }

    long long DynamicMST::mst_weight() {
        long long s = 0;
        for(int u = 0; u < n; ++u)
            if(nodes[u].p.v >= 0) s += nodes[u].p.w;
        return s;
    }

    void DynamicMST::print_mst() {
        for(int i = 1556; i < n; ++i)
            if (i == 1556 || i == 1732) {
                printf("%d: p=(%d,%d,%d), min_e=(%d,%d,%d), sub_cnt=%d\n", i, nodes[i].p.u, nodes[i].p.v,
                       nodes[i].p.w, nodes[i].min_e.u, nodes[i].min_e.v, nodes[i].min_e.w, nodes[i].sub_cnt);
                if (i == 1732 && nodes[i].p.v != 1556 && nodes[1556].p.v != 1732) {
                    exit(0);
                }
            }
        printf( "=====\n" );
    }

    bool DynamicMST::is_decedent(int u, int v) {
        for(; u != -1; u = nodes[u].p.v)
            if(u == v) return true;
        return false;
    }

    bool DynamicMST::is_out_edge(int u, Edge e) {
        return is_decedent(e.u, u) && !is_decedent(e.v,u);
    }


    bool DynamicMST::get_temporal_edge(char *line, int &a, int &b, int &t, int num_cnt) {
        if( !isdigit(line[0]) ) return false;
        vector<char*> v_num;
        int len = (int) strlen(line);
        for( int i = 0; i < len; ++i )
            if( !isdigit(line[i]) && line[i] != '.') line[i] = '\0';
            else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
        if( (int) v_num.size() != num_cnt ) return false;
        sscanf( v_num[0], "%d", &a );
        sscanf( v_num[1], "%d", &b );
        sscanf( v_num[num_cnt-1], "%d", &t );
        return true;
    }

    void DynamicMST::create_stream(string path, bool is_bipartite) {
        FILE *fin = fopen((path + "graph.txt").c_str(), "r" );
        char line[MAXST];
        int n = 0, n1 = 0, n2 = 0, a, b, num_cnt = get_num_cnt(path + "graph.txt"), t;
        vector<pair<pair<int,int>, pair<int,int> > > el;
        long long m = 0;

        printf("Loading text, num_cnt = %d, is_bipartite = %s...\n", num_cnt, is_bipartite?"true":"false");
        while(fgets( line, MAXST, fin)) {
            if(!get_temporal_edge(line, a, b, t, num_cnt)) continue;
            if(a < 0 || b < 0 || ((!is_bipartite) && a == b)) continue;
            el.push_back(make_pair(make_pair(t,rand()%1000000),make_pair(a, b)));
            n = max(max(n, a+1), b+1);
            n1 = max(n1, a+1);
            n2 = max(n2, b+1);
            if((++m) % (long long) 10000000 == 0) printf("%lld lines finished\n", m);
        }
        fclose(fin);

        if(is_bipartite) {
            printf( "n1 = %d, n2 = %d\n", n1, n2 );
            n = n1 + n2;
            for(int i = 0; i < (int) el.size(); ++i) el[i].second.second += n1;
        }

        for(int i = 0; i < (int) el.size(); ++i)
            if(el[i].second.first > el[i].second.second)
                swap(el[i].second.first, el[i].second.second);

        sort(el.begin(),el.end());

        m = 0;
        for(int i = 0; i < (int) el.size(); ++i)
            if(i == 0 || el[i] != el[i-1])
                el[m++] = el[i];

        FILE *fout = fopen((path + "graph.stream").c_str(), "wb");
        fwrite(&n, sizeof(int), 1, fout);
        fwrite(&m, sizeof(long long), 1, fout);
        fwrite(el.data(), sizeof(pair<pair<int,int>,pair<int,int> >), m, fout);

        fclose(fout);
        printf("Created stream file, n = %d, m = %lld\n", n, m);
    }
}