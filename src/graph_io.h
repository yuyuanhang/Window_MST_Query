#ifndef GRAPH_IO_H
#define GRAPH_IO_H

#include "data_type.h"
#include <string>
#include <vector>

class GraphIO {
public:
    GraphIO(const std::string& filename) : filename_(filename) {}

    std::vector<Edge> read_graph();
    void write_graph(const std::string &filename, std::vector<Edge> edges);

private:
    std::string filename_;
};

#endif // GRAPH_IO_H