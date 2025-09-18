#include "graph_io.h"
#include "graph.h"
#include "mst.h"
#include "window_generator.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>

enum class CommandType {
    GenerateQuery,
    BuildIndex,
    BuildIndexPara,
    BuildIndexInc,
    ExecuteQuery,
    ValidateIndex,
    QTCount,
    Unknown
};

std::map<std::string, CommandType> commandMap = {
        {"--generate-query", CommandType::GenerateQuery},
        {"-gq", CommandType::GenerateQuery},
        {"--build-index", CommandType::BuildIndex},
        {"-bi", CommandType::BuildIndex},
        {"--build-index-parallel", CommandType::BuildIndexPara},
        {"-bip", CommandType::BuildIndexPara},
        {"--build-index-incremental", CommandType::BuildIndexInc},
        {"-bii", CommandType::BuildIndexInc},
        {"--execute-query", CommandType::ExecuteQuery},
        {"-eq", CommandType::ExecuteQuery}
};

CommandType parseCommand(const std::string& cmd) {
    if (commandMap.find(cmd) != commandMap.end()) {
        return commandMap[cmd];
    }
    return CommandType::Unknown;
}

void help() {
    std::cout << "Usage: [command] [options]" << std::endl;
    std::cout << "--generate-query/-gq [data graph] [window_ratio] [num_queries] [query]" << std::endl;
    std::cout << "--build-index/-bi [data graph] [method_type] [index]" << std::endl;
    std::cout << "--execute-query/-eq [data graph] [query] [method_type] [index]" << std::endl;
    std::cout << "--build-index-parallel/-bip [data graph] [num_thread] [index]" << std::endl;
    std::cout << "--build-index-incremental/-bii [data graph] [edge_ratio]" << std::endl;
}

void generate_queries(const std::string &tgraph, const float window_ratio, const int num_queries, std::string query) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges);
    // Modify the graph data
    graph.group_edges_by_timestamp();
    TimeWindowGenerator window_generator(graph.get_num_timestamps(), query);
    window_generator.generate(num_queries, window_ratio);
    window_generator.save_windows();
}

void build_index(const std::string &tgraph, const std::string &type, const std::string &index) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges);
    // Modify the graph data
    graph.group_edges_by_timestamp();
    MSTIndex mst_index(&graph);

    if (type == "base")
    {
        std::string bl_index(index+"_bl.bin");
        mst_index.create_bl_idx_plus();
        mst_index.save_bl_idx(bl_index);
    } else if (type == "op")
    {
        std::string bl_index(index+"_bl.bin");
        mst_index.create_bl_idx_pro();
        mst_index.save_bl_idx(bl_index);
    } else if (type == "d")
    {
        std::string d_index(index+"_d.bin");
        mst_index.create_bl_idx_pro();
        mst_index.create_qt_idx();
        mst_index.save_qt_idx(d_index);
    } else if (type == "b")
    {
        std::string b_index(index+"_b.bin");
        mst_index.create_bl_idx_pro();
        mst_index.create_hy_idx(true);
        mst_index.save_hy_idx(b_index);
    }
}

void build_index_para(const std::string &tgraph, const int num_threads, const std::string &index) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges);
    // Modify the graph data
    graph.group_edges_by_timestamp();
    MSTIndex mst_index(&graph);
    std::string bl_index(index + "_bl_para_" + std::to_string(num_threads) + ".bin");
    mst_index.create_bl_idx_pro_parallel(num_threads);
    mst_index.save_bl_idx(bl_index);
}

void build_index_inc(const std::string &tgraph, const float edge_ratio) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges, edge_ratio);
    // Modify the graph data
    graph.group_edges_by_timestamp();
    MSTIndex mst_index(&graph);
    mst_index.time_window_sz = graph.nb_tbase;
    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds delete_time(0);
    mst_index.create_bl_idx_pro();
    for (auto i = 0; i < graph.nb_tstream; i++)
    {
        mst_index.update_bl_idx();
        mst_index.create_hy_idx(true);

        auto d0 = std::chrono::high_resolution_clock::now();
        mst_index.delete_hy_idx();
        auto d1 = std::chrono::high_resolution_clock::now();
        delete_time += std::chrono::duration_cast<std::chrono::microseconds>(d1 - d0);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::chrono::microseconds final_time = dur - delete_time;
    std::cout << "bl index construction time (pro && inc): "
              << (double)final_time.count() / 1000 << " ms" << std::endl;
}

void execute_queries(const std::string &tgraph, std::string query, const std::string &type, const std::string &index) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges);
    // Modify the graph data
    graph.group_edges_by_timestamp();

    TimeWindowGenerator window_generator(graph.get_num_timestamps(), query);
    window_generator.load_windows();

    MSTIndex mst_index(&graph);
    if (type == "online")
    {
        mst_index.query_online_global(window_generator.windows_);
    } else if (type == "base")
    {
        std::string bl_index(index+"_bl.bin");
        mst_index.load_bl_idx(bl_index);
        mst_index.query_bl(window_generator.windows_);
    } else if (type == "d")
    {
        std::string d_index(index+"_d.bin");
        mst_index.load_qt_idx(d_index);
        mst_index.query_qt(window_generator.windows_);
    } else if (type == "b")
    {
        std::string b_index(index+"_b.bin");
        mst_index.load_hy_idx(b_index);
        mst_index.query_hy(window_generator.windows_);
    }
}

int main(int argc, char* argv[])
{
    if (argc < 2) { help(); return 0; }

    CommandType cmdType = parseCommand(argv[1]);
    switch (cmdType) {
        case CommandType::GenerateQuery: {
            std::cout << "Generating query..." << std::endl;
            if (argc < 6) { help(); return 0; }
            std::string tgraph = std::string(argv[2]);
            float window_ratio = std::stof(argv[3]);
            int num_queries = std::stoi(argv[4]);
            std::string query = std::string(argv[5]);

            generate_queries(tgraph, window_ratio, num_queries, query);
            break;
        }
        case CommandType::BuildIndex: {
            std::cout << "Building index..." << std::endl;
            if (argc < 5) { help(); return 0; }
            std::string tgraph = std::string(argv[2]);
            std::string type = std::string(argv[3]);
            std::string index = std::string(argv[4]);
            build_index(tgraph, type, index);
            break;
        }
        case CommandType::BuildIndexPara: {
            std::cout << "Building index in parallel..." << std::endl;
            if (argc < 5) { help(); return 0; }
            std::string tgraph = std::string(argv[2]);
            int num_threads = std::stoi(argv[3]);
            std::string index = std::string(argv[4]);
            build_index_para(tgraph, num_threads, index);
            break;
        }
        case CommandType::BuildIndexInc: {
            std::cout << "Updating index..." << std::endl;
            if (argc < 4) { help(); return 0; }
            std::string tgraph = std::string(argv[2]);
            float edge_ratio = std::stof(argv[3]);
            build_index_inc(tgraph, edge_ratio);
            break;
        }
        case CommandType::ExecuteQuery: {
            std::cout << "Executing query..." << std::endl;
            if (argc < 6) { help(); return 0; }
            std::string tgraph = std::string(argv[2]);
            std::string query = std::string(argv[3]);
            std::string type = std::string(argv[4]);
            std::string index = std::string(argv[5]);

            execute_queries(tgraph, query, type, index);
            break;
        }
        default: {
            std::cout << "Unknown command." << std::endl;
            help();
            break;
        }
    }
}
