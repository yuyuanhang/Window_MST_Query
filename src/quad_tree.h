//
// Created by yuanyu on 12/06/24.
//

#ifndef QUAD_TREE_H
#define QUAD_TREE_H

#include "data_type.h"
#include <iostream>
#include <vector>
#include <set>
#include <memory>
#include <climits>
#include <fstream>

// x: t1 <= te < t2; y: s < ts <= t
struct Rectangle {
    int x, y, width, height;

    Rectangle() : x(0), y(0), width(0), height(0) {}

    Rectangle(int x, int y, int width, int height)
            : x(x), y(y), width(width), height(height) {}

    bool intersects(const Rectangle& other) const {
        return !(x > other.x + other.width || x + width < other.x ||
                 y > other.y + other.height || y + height < other.y);
    }

    Rectangle intersection(const Rectangle& other) const {
        int newX = std::max(x, other.x);
        int newY = std::max(y, other.y);
        int newWidth = std::min(x + width, other.x + other.width) - newX;
        int newHeight = std::min(y + height, other.y + other.height) - newY;
        if (newWidth > 0 && newHeight > 0) {
            return Rectangle(newX, newY, newWidth, newHeight);
        }
        return Rectangle(0, 0, 0, 0);
    }

    bool isEmpty() const {
        return width <= 0 || height <= 0;
    }

    bool operator==(const Rectangle& other) const {
        return x == other.x && y == other.y && width == other.width && height == other.height;
    }
};

struct Point {
    int x, y;

    Point(int x, int y)
            : x(x), y(y) {}

    bool contained(const Rectangle& other) const {
        return !(x >= other.x + other.width || x < other.x ||
                 y > other.y + other.height || y <= other.y);
    }
};

static const int LEFT_TOP = 0;
static const int LEFT_BOTTOM = 1;
static const int RIGHT_TOP = 2;
static const int RIGHT_BOTTOM = 3;

class QuadTree {
public:
    QuadTree(const Rectangle& boundary) {
        root = std::make_unique<Node>(boundary);
    }

    bool insert(const Rectangle& t, Edge e) {
        return root->insert(std::make_pair(t, e));
    }

    int count(const Rectangle& r) {
        return root->count(r);
    }

    void retrieve(const Point& p, std::vector<Edge> &res) {
        root->retrieve(p, res);
    }

    void serialize(const std::string& filename) const {
        std::ofstream ofs(filename, std::ios::binary);
        if (ofs.is_open()) {
            root->serialize(ofs);
            ofs.close();
        }
    }

    void deserialize(const std::string& filename) {
        std::ifstream ifs(filename, std::ios::binary);
        if (ifs.is_open()) {
            root = std::make_unique<Node>(Rectangle(0, 0, 0, 0)); // root node
            root->deserialize(ifs);
            ifs.close();
        }
    }

    void clear() {
        root->clear();
    }

private:
    struct Node {
        Rectangle boundary;
        std::vector<Edge> bound_edges;
        std::vector<std::unique_ptr<Node>> childs;
        bool childEmpty[4];

        Node(const Rectangle& boundary)
                : boundary(boundary) {}

        bool insert(std::pair<Rectangle, Edge> t) {
            auto sub_t = t.first.intersection(boundary);
            // non-overlap
            if (sub_t.isEmpty()) {
                return false;
            }
            // node is covered by the rectangle
            if (sub_t == boundary) {
                // insert this rectangle into this node
                bound_edges.push_back(std::move(t.second));
            } else { // some child nodes are covered by the rectangle
                if (childs.empty()) {
                    split();
                }
                auto quadIndex = quadrantOf(t.first);
                for (auto &childIndex : quadIndex) {
                    childs[childIndex]->insert(std::move(t));
                }
            }

            return true;
        }

        int count(Rectangle r) {
            auto sub_r = r.intersection(boundary);
            // non-overlap
            if (sub_r.isEmpty()) {
                return 0;
            }
            // node is covered by the rectangle
            if (sub_r == boundary) {
                return 1;
            } else { // some child nodes are covered by the rectangle
                if (childs.empty()) {
                    split();
                    for (auto i = 0; i < 4; i++) {
                        if (childs[i]->boundary.isEmpty()) {
                            childEmpty[i] = true;
                        } else {
                            childEmpty[i] = false;
                        }
                    }
                }
                int sum = 0;
                auto quadIndex = quadrantOf(r);
                for (auto &childIndex : quadIndex) {
                    if (!childEmpty[childIndex]) {
                        sum += childs[childIndex]->count(r);
                    }
                }
                return sum;
            }
        }

        void split() {
            if (childs.empty()) {
                int hmid = boundary.x + boundary.width / 2;
                int vmid = boundary.y + boundary.height / 2;
        
                int halfWidth1 = boundary.width / 2;
                int halfHeight1 = boundary.height / 2;
                int halfWidth2 = boundary.width - halfWidth1;
                int halfHeight2 = boundary.height - halfHeight1;
        
                childs.resize(4);
                childs[LEFT_TOP] = std::make_unique<Node>(Rectangle(boundary.x, boundary.y, halfWidth1, halfHeight1));
                childs[LEFT_BOTTOM] = std::make_unique<Node>(Rectangle(boundary.x, vmid, halfWidth1, halfHeight2));
                childs[RIGHT_TOP] = std::make_unique<Node>(Rectangle(hmid, boundary.y, halfWidth2, halfHeight1));
                childs[RIGHT_BOTTOM] = std::make_unique<Node>(Rectangle(hmid, vmid, halfWidth2, halfHeight2));
            }
        }

        [[nodiscard]] std::vector<int> quadrantOf(const Rectangle& rect) const {
            std::vector<int> res;

            if (rect.intersects(childs[LEFT_TOP]->boundary)) {
                res.push_back(LEFT_TOP);
            } if (rect.intersects(childs[LEFT_BOTTOM]->boundary)) {
                res.push_back(LEFT_BOTTOM);
            } if (rect.intersects(childs[RIGHT_TOP]->boundary)) {
                res.push_back(RIGHT_TOP);
            } if (rect.intersects(childs[RIGHT_BOTTOM]->boundary)) {
                res.push_back(RIGHT_BOTTOM);
            }
            return res;
        }

        void retrieve(const Point& p, std::vector<Edge> &resultSet) const {
            resultSet.insert(resultSet.end(), bound_edges.begin(), bound_edges.end());
            if (!childs.empty()) {
                for (auto i = 0; i < 4; i++) {
                    if (p.contained(childs[i]->boundary)) {
                        childs[i]->retrieve(p, resultSet);
                    }
                }
            }
        }

        void serialize(std::ofstream& ofs) const {
            ofs.write(reinterpret_cast<const char*>(&boundary), sizeof(boundary));

            int edgeCount = bound_edges.size();
            ofs.write(reinterpret_cast<const char*>(&edgeCount), sizeof(edgeCount));
            for (auto& edge : bound_edges) {
                ofs.write(reinterpret_cast<const char*>(&edge), sizeof(Edge));
            }
            int hasChildren = !childs.empty();
            ofs.write(reinterpret_cast<const char*>(&hasChildren), sizeof(hasChildren));
            if (hasChildren) {
                for (const auto& child : childs) {
                    child->serialize(ofs);
                }
            }
        }

        void deserialize(std::ifstream& ifs) {
            ifs.read(reinterpret_cast<char*>(&boundary), sizeof(boundary));

            int edgeCount;
            ifs.read(reinterpret_cast<char*>(&edgeCount), sizeof(edgeCount));
            bound_edges.resize(edgeCount);
            for (int i = 0; i < edgeCount; ++i) {
                ifs.read(reinterpret_cast<char*>(&bound_edges[i]), sizeof(Edge));
            }
            int hasChildren;
            ifs.read(reinterpret_cast<char*>(&hasChildren), sizeof(hasChildren));
            if (hasChildren) {
                childs.resize(4);
                for (int i = 0; i < 4; ++i) {
                    childs[i] = std::make_unique<Node>(Rectangle(0, 0, 0, 0));
                    childs[i]->deserialize(ifs);
                }
            }
        }

        void clear() {
            bound_edges.clear();
            bound_edges.shrink_to_fit();
            if (!childs.empty()) {
                for (auto i = 0; i < 4; i++) {
                    childs[i]->clear();
                }
            }
        }
    };

    std::unique_ptr<Node> root;
};

#endif //QUAD_TREE_H
