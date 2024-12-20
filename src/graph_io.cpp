#include "graph_io.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<Edge> GraphIO::read_graph() {
    std::vector<Edge> edges;
    std::ifstream infile(filename_, std::ios::binary);
    if (infile.is_open()) {
        std::cout << "File opened successfully." << std::endl;
    } else {
        std::cout << "Failed to open file." << std::endl;
    }
    int num_edges = 0;
    infile.read(reinterpret_cast<char*>(&num_edges), sizeof(int));
    edges.resize(num_edges);
    infile.read(reinterpret_cast<char*>(edges.data()), num_edges * sizeof(Edge));
    infile.close();
    /*std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int src_id, dest_id;
        bint timestamp;
        float weight = 0;

        if (!(iss >> src_id >> dest_id >> weight >> timestamp)) {
            throw std::runtime_error("Error parsing line: " + line);
        }

        edges.emplace_back(src_id, dest_id, weight, timestamp);
    }*/
    return edges;
}

void GraphIO::write_graph(const std::string &filename, std::vector<Edge> edges) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (auto edge : edges) {
            outFile << edge.src_id << " " << edge.dest_id << " " << edge.weight << " " << edge.timestamp << "\n";
        }
        outFile.close();
        std::cout << "Edges have been written to file " << filename << std::endl;
    } else {
        std::cout << "Unable to open file " << filename << std::endl;
    }
}