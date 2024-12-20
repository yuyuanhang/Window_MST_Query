
#include "DynamicMST.h"

int get_dep(DynamicMST &g, int u, vector<int> &dep) {
	if(u == -1) return -1;
	if(dep[u] != -1) return dep[u];
	return (dep[u] = get_dep(g, g.nodes[u].p.v, dep) + 1);
}

void print_dep(DynamicMST &g) {
	vector<int> dep(g.n,-1);
	int max_dep = 0;
	long long sum_dep = 0;
	for(int u = 0; u < g.n; ++u) {
		sum_dep += get_dep(g, u, dep);
		max_dep = max(max_dep, dep[u]);
	}
	printf( "max_dep = %d, avg_dep = %0.3lf\n", max_dep, 1.0*sum_dep/g.n);
}

void preprocess_unweighted(string path) {
	printf( "Preprocessing unweighted graph %s\n", path.c_str() );
	FILE *fin = fopen((path+"graph-org.txt").c_str(), "r");
	FILE *fout = fopen((path+"graph.txt").c_str(), "w");

	int a, b, num_cnt = DynamicMST::get_num_cnt(path+"graph-org.txt");

	char line[MAXST];

	vector<pair<int,int> > e_all;

	while(fgets( line, MAXST, fin )) {
		if( !isdigit(line[0]) ) continue;
		vector<char*> v_num;
		int len = (int) strlen(line);
		for( int i = 0; i < len; ++i )
			if( !isdigit(line[i]) && line[i] != '.') line[i] = '\0';
			else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
		if((int) v_num.size() != num_cnt) continue;
		sscanf(v_num[0], "%d", &a);
		sscanf(v_num[1], "%d", &b);
		if(a > b) swap(a, b);
		if(a != b) e_all.push_back(make_pair(a,b));
	}

	sort(e_all.begin(), e_all.end());
	long long m = 0;

	for(size_t i = 0; i < e_all.size(); ++i)
		if(i == 0 || e_all[i] != e_all[i-1]) {
			fprintf(fout, "%d %d %d\n", e_all[i].first, e_all[i].second, rand()/64+1 );
			++m;
		}

	printf( "m=%lld\n", m );

	fclose(fin);
	fclose(fout);
}

void sample(string path, long long n_edges) {
	DynamicMST g(path);
	clock_t t = clock();
	g.init();
	printf( "Initialize time = %0.3lf sec\n", (clock()-t)*1.0/CLOCKS_PER_SEC);

	vector<Edge> l;
	g.sample_edges(l, n_edges);
	printf( "Sampled edges = %lld\n", (long long)l.size());
	print_dep(g);
	t = clock();
	long long c0 = 0, c1 = 0, c2 = 0;
	for(long long i = 0; i < (long long)l.size(); ++i) {

		int r = g.delete_edge(l[i].u,l[i].v,l[i].w);
		if(r == 0) ++c0;
		if(r == 1) ++c1;
		if(r == 2) ++c2;
	}

	printf( "Deletion time = %0.3lf sec, c0 = %lld, c1 = %lld, c2 = %lld\n", (clock()-t)*1.0/CLOCKS_PER_SEC, c0, c1, c2);
	printf( "s1 = %lld, s2 = %lld, s3 = %lld, t = %lld, t3 = %lld, avg_s1 = %0.3lf, avg_s2 = %0.3lf, avg_s3 = %0.3lf\n",
			g.s1, g.s2, g.s3, g.t, g.t3, g.s1*1.0/g.t, g.s2*1.0/g.t, g.s3*1.0/g.t3);
	printf( "MST weight = %lld\n", g.mst_weight() );

	t = clock();
	c0 = 0; c1 = 0; c2 = 0;

	for(long long i = 0; i < (long long)l.size(); ++i) {
		int r = g.insert_edge(l[i].u,l[i].v,l[i].w);
		if(r == 0) ++c0;
		if(r == 1) ++c1;
		if(r == 2) ++c2;
	}

	printf( "Insertion time = %0.3lf sec, c0 = %lld, c1 = %lld, c2 = %lld\n", (clock()-t)*1.0/CLOCKS_PER_SEC, c0, c1, c2);
	printf( "MST weight = %lld\n", g.mst_weight() );
	printf("cand = %0.3lf, qq = %0.3lf, nontree = %0.3lf\n", candsize * 1.0 / candtimes, qqsize * 1.0 / qqtimes, nontreesize * 1.0 / nontreetimes);
	printf("einf = %0.3lf, qcheck = %0.3lf\n", einf * 1.0 / qqtimes, qcheck * 1.0 / qqtimes);
	printf("aa = %0.3lf, bb = %0.3lf, cc = %0.3lf, dd = %0.3lf\n", aa * 1.0 / qqtimes, bb * 1.0 / qqtimes, cc * 1.0 / qqtimes, dd * 1.0 / qqtimes);
	print_dep(g);
}

void stream(string path, int window_size) {
	clock_t t = clock();
	DynamicMST g(path, false);
	g.init();
	printf( "Initialize time = %0.3lf sec\n", (clock()-t)*1.0/CLOCKS_PER_SEC);

	vector<pair<pair<int,int>,pair<int,int> > > l;
	int n;
	long long m;

	FILE *fin = fopen( (path+"graph.stream").c_str(), "rb" );
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(long long), 1, fin);
	printf( "n = %d, m = %lld\n", n, m );

	l.resize(m);
	fread(l.data(),sizeof(pair<pair<int,int>,pair<int,int> >), m, fin);

	fclose(fin);

	int ws = window_size;
	if(window_size>0 && window_size<100)
	{
		int tmin = INT_MAX;
		int tmax = 0;
		for (long long j = 0; j < (long long)l.size(); ++j)
		{
			tmin = min(tmin, l[j].first.first);
			tmax = max(tmax, l[j].first.first);
		}
		double w = double(window_size)/100;
		window_size = (int)((tmax - tmin) * w);
		cout << "tmin = " << tmin << " tmax = " << tmax << " window_size = " << window_size <<" "<<w<< endl;
	}
	if(window_size < 0) window_size = l[l.size()-1].first.first;

	t = clock();
	int c0_ins = 0, c1_ins = 0, c2_ins = 0, c0_del = 0, c1_del = 0, c2_del = 0;
	for(int s = 0, i = 0; i < (int) l.size(); ++i) {
		int r = g.insert_edge(l[i].second.first, l[i].second.second, l[i].first.second);
		if(r == 0) ++c0_ins;
		if(r == 1) ++c1_ins;
		if(r == 2) ++c2_ins;

		if(i % 1000000 == 0) printf( "%d/%d\n", i, (int) l.size());

		for(; l[i].first.first-l[s].first.first > window_size; ++s) {
			r = g.delete_edge(l[s].second.first, l[s].second.second, l[s].first.second);
			if(r == 0) ++c0_del;
			if(r == 1) ++c1_del;
			if(r == 2) ++c2_del;
		}
	}


	printf( "Update time = %0.3lf sec, c0_ins = %d, c1_ins = %d, c2_ins = %d, c0_del = %d, c1_del = %d, c2_del = %d, n_ins = %d, n_del = %d\n",
			(clock()-t)*1.0/CLOCKS_PER_SEC, c0_ins, c1_ins, c2_ins, c0_del, c1_del, c2_del, c0_ins+c1_ins+c2_ins, c0_del+c1_del+c2_del);
	print_dep(g);
	if(ws == -2) {
		t = clock();
		for(int s = 0, i = 0; i < (int) l.size(); ++i) {
			int r = g.delete_edge(l[i].second.first, l[i].second.second, l[i].first.second);
			if(r == 0) ++c0_del;
			if(r == 1) ++c1_del;
			if(r == 2) ++c2_del;
		}
		printf( "Deletion time = %0.3lf sec, c0_del = %d, c1_del = %d, c2_del = %d, n_del = %d\n",
					(clock()-t)*1.0/CLOCKS_PER_SEC, c0_del, c1_del, c2_del, c0_del+c1_del+c2_del);
	}

	printf( "s1 = %lld, s2 = %lld, s3 = %lld, t = %lld, t3 = %lld, avg_s1 = %0.3lf, avg_s2 = %0.3lf, avg_s3 = %0.3lf\n",
			g.s1, g.s2, g.s3, g.t, g.t3, g.s1*1.0/g.t, g.s2*1.0/g.t, g.s3*1.0/g.t3);
	printf("cand = %0.3lf, qq = %0.3lf, nontree = %0.3lf\n", candsize * 1.0 / candtimes, qqsize * 1.0 / qqtimes, nontreesize * 1.0 / nontreetimes);
	printf("einf = %0.3lf, qcheck = %0.3lf\n", einf * 1.0 / qqtimes, qcheck * 1.0 / qqtimes);
	printf("aa = %0.3lf, bb = %0.3lf, cc = %0.3lf, dd = %0.3lf\n", aa * 1.0 / qqtimes, bb * 1.0 / qqtimes, cc * 1.0 / qqtimes, dd * 1.0 / qqtimes);
}


void preprocess_road(string path) {
	printf( "Preprocessing road %s\n", path.c_str() );
	FILE *fin = fopen((path+"graph-org.txt").c_str(), "r");
	FILE *fout = fopen((path+"graph.txt").c_str(), "w");
	char line[MAXST];

	while(fgets( line, MAXST, fin )) {
		if(line[0] != 'a') continue;
		vector<char*> v_num;
		int len = (int) strlen(line);
		for(int i = 0; i < len; ++i)
			if(!isdigit(line[i]) && line[i] != '.') line[i] = '\0';
			else if(i>0 && !line[i-1]) v_num.push_back(line+i);
		if(v_num.size() != 3) printf( "#" );
		int u = atoi(v_num[0]), v = atoi(v_num[1]), w = atoi(v_num[2]);
		if(u<v) fprintf( fout, "%d %d %d\n", u, v, w);

	}
	fclose(fin);
	fclose(fout);
}

int main(int argc, char *argv[]) {
	printf( "argc=%d\n", argc );
	for( int i = 0; i < argc; ++i )
		printf( "argv[%d]=%s\n", i, argv[i] );

	setvbuf(stdout, NULL, _IONBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);
	clock_t t = clock();

	if( argc > 1 ) {
		if(strcmp(argv[1], "txt-to-bin") == 0) {
			bool multi_edge = false;
			if(argc > 3) multi_edge = atoi(argv[3]);
			DynamicMST::create_bin( /*dataset*/argv[2], multi_edge );
		} else if( strcmp( argv[1], "sample" ) == 0 )
			sample(argv[2], atoll(argv[3]));
		else if( strcmp( argv[1], "preprocess-unweighted" ) == 0 )
			preprocess_unweighted(argv[2]);
		else if( strcmp(argv[1], "txt-to-stream" ) == 0 ) {
			bool is_bipartite = false;
			if(argc > 3 && strcmp(argv[3], "bipartite") == 0) is_bipartite = true;
			DynamicMST::create_stream( /*dataset*/argv[2], is_bipartite );
		} else if( strcmp( argv[1], "stream" ) == 0 )
			stream(argv[2], atoi(argv[3]));
		else if( strcmp( argv[1], "preprocess-road" ) == 0 )
			preprocess_road(argv[2]);
	}

	if( argc <= 1) {
	}

	t = clock() - t;
	printf( "Total time=%0.3lf seconds\n", t*1.0/CLOCKS_PER_SEC);

	return 0;
}
