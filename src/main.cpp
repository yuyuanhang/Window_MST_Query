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
        {"--execute-query", CommandType::ExecuteQuery},
        {"-eq", CommandType::ExecuteQuery},
        {"--validate-index", CommandType::ValidateIndex},
        {"-vi", CommandType::ValidateIndex},
        {"--qt-count", CommandType::QTCount},
        {"-qc", CommandType::QTCount}
};

CommandType parseCommand(const std::string& cmd) {
    if (commandMap.find(cmd) != commandMap.end()) {
        return commandMap[cmd];
    }
    return CommandType::Unknown;
}

void help() {
    std::cout << "Usage: [command] [options]" << std::endl;
    std::cout << "--generate-query/-gq [output] [num_queries] [data graph] [window_ratio]" << std::endl;
    std::cout << "--build-index/-bi [output] [data graph] [edge_ratio]" << std::endl;
    std::cout << "--execute-query/-eq [output] [queries] [data graph] [index] [method_type]" << std::endl;
    std::cout << "--validate-index/-vi [queries] [data graph] [index] [window_ratio]" << std::endl;
    std::cout << "--qt-count/-qc [output] [data graph] [edge_ratio]" << std::endl;
}

void generate_queries(int num_queries, std::string query, const std::string &tgraph, float window_ratio) {
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

void build_index(const std::string &tgraph, std::string index, float edge_ratio) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges, edge_ratio);
    // Modify the graph data
    graph.group_edges_by_timestamp();
    MSTIndex mst_index(&graph);
    //mst_index.create_bl_idx();
    //mst_index.save_bl_idx(index);
    std::string ratio;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << edge_ratio;
    ratio = oss.str();
    std::string bl_index(index+"_bl_"+ratio+".bin");
    mst_index.create_bl_idx_plus();
    /*mst_index.save_bl_idx(bl_index);

    std::string op_index(index+"_op_"+ratio+".bin");
    mst_index.create_op_idx();
    mst_index.save_op_idx(op_index);

    std::string hy1_index(index+"_hy1_"+ratio+".bin");
    mst_index.create_hy_idx(true);
    mst_index.save_hy_idx(hy1_index);

    std::string qt_index(index+"_qt_"+ratio+".bin");
    mst_index.create_qt_idx();
    mst_index.save_qt_idx(qt_index);*/

    /*std::string hy0_index(index+"_hy0_"+ratio+".bin");
    mst_index.create_hy_idx(false);
    mst_index.save_hy_idx(hy0_index);*/
}

void qt_count(const std::string &tgraph, std::string index, float edge_ratio) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges, edge_ratio);
    // Modify the graph data
    graph.group_edges_by_timestamp();
    MSTIndex mst_index(&graph);
    mst_index.qt_idx_count(index);
}

void execute_queries(std::string query, const std::string &tgraph, std::string index, std::string ratio, const std::string &res) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges);
    // Modify the graph data
    graph.group_edges_by_timestamp();

    TimeWindowGenerator window_generator(graph.get_num_timestamps(), query);
    window_generator.load_windows();

    MSTIndex mst_index(&graph);
    /*mst_index.query_online(window_generator.windows_);
    std::string bl_index(index+"_bl_"+ratio+".bin");
    mst_index.load_bl_idx(bl_index);
    mst_index.query_bl(window_generator.windows_);*/

    std::string op_index(index+"_op_"+ratio+".bin");
    mst_index.load_op_idx(op_index);
    mst_index.query_op(window_generator.windows_);
    
    /*std::string qt_index(index+"_qt_"+ratio+".bin");
    mst_index.load_qt_idx(qt_index);
    mst_index.query_qt(window_generator.windows_);*/
    
    std::string hy1_index(index+"_hy1_"+ratio+".bin");
    mst_index.load_hy_idx(hy1_index);
    mst_index.query_hy(window_generator.windows_);
}

void validate_index(std::string query, const std::string &tgraph, std::string index, std::string ratio) {
    GraphIO graph_io(tgraph);
    // Read the graph data from file
    std::vector<Edge> edges = graph_io.read_graph();
    TemporalGraph graph(edges);
    // Modify the graph data
    graph.group_edges_by_timestamp();

    TimeWindowGenerator window_generator(graph.get_num_timestamps(), query);
    window_generator.load_windows();

    MSTIndex mst_index(&graph);
    
    mst_index.validate(window_generator.windows_, index, ratio);
}

int main(int argc, char* argv[])
{
    if (argc < 2) { help(); return 0; }

    CommandType cmdType = parseCommand(argv[1]);
    switch (cmdType) {
        case CommandType::GenerateQuery: {
            std::cout << "Generating query..." << std::endl;
            if (argc < 6) { help(); return 0; }
            std::string query = std::string(argv[2]);
            int num_queries;
            try {
                num_queries = std::stoi(argv[3]);
            } catch (const std::invalid_argument &e) {
                std::cerr << "Invalid number: " << argv[3] << std::endl;
            } catch (const std::out_of_range &e) {
                std::cerr << "Number out of range: " << argv[3] << std::endl;
            }
            std::string tgraph = std::string(argv[4]);
            float window_ratio;
            try {
                window_ratio = std::stof(argv[5]);
            } catch (const std::invalid_argument &e) {
                std::cerr << "Invalid number: " << argv[5] << std::endl;
            } catch (const std::out_of_range &e) {
                std::cerr << "Number out of range: " << argv[5] << std::endl;
            }
            generate_queries(num_queries, query, tgraph, window_ratio);
            break;
        }
        case CommandType::BuildIndex: {
            std::cout << "Building index..." << std::endl;
            if (argc < 5) { help(); return 0; }
            std::string index = std::string(argv[2]);
            std::string tgraph = std::string(argv[3]);
            float edge_ratio;
            try {
                edge_ratio = std::stof(argv[4]);
            } catch (const std::invalid_argument &e) {
                std::cerr << "Invalid number: " << argv[4] << std::endl;
            } catch (const std::out_of_range &e) {
                std::cerr << "Number out of range: " << argv[4] << std::endl;
            }
            build_index(tgraph, index, edge_ratio);
            break;
        }
        case CommandType::ExecuteQuery: {
            std::cout << "Executing query..." << std::endl;
            if (argc < 5) {
                help();
                return 0;
            }
            std::string res = std::string(argv[2]);
            std::string query = std::string(argv[3]);
            std::string tgraph = std::string(argv[4]);
            std::string index;
            std::string ratio;
            if (argc == 7) {
                index = std::string(argv[5]);
                ratio = std::string(argv[6]);
            }
            execute_queries(query, tgraph, index, ratio, res);
            break;
        }
        case CommandType::ValidateIndex: {
            std::cout << "Validating index..." << std::endl;
            if (argc < 6) {
                help();
                return 0;
            }
            std::string query = std::string(argv[2]);
            std::string tgraph = std::string(argv[3]);
            std::string index = std::string(argv[4]);
            std::string edge_ratio = std::string(argv[5]);
            validate_index(query, tgraph, index, edge_ratio);
            break;
        }
        case CommandType::QTCount: {
            std::cout << "counting quad-tree statistics..." << std::endl;
            if (argc < 5) {
                help();
                return 0;
            }
            std::string index = std::string(argv[2]);
            std::string tgraph = std::string(argv[3]);
            float edge_ratio;
            try {
                edge_ratio = std::stof(argv[4]);
            } catch (const std::invalid_argument &e) {
                std::cerr << "Invalid number: " << argv[4] << std::endl;
            } catch (const std::out_of_range &e) {
                std::cerr << "Number out of range: " << argv[4] << std::endl;
            }
            qt_count(tgraph, index, edge_ratio);
            break;
        }
        default: {
            std::cout << "Unknown command." << std::endl;
            help();
            break;
        }
    }
}
