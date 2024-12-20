//
// Created by yuanyu on 11/07/24.
//

#ifndef HYBRID_H
#define HYBRID_H

#include "data_type.h"
#include "quad_tree.h"
#include "mst.h"
#include <iostream>
#include <memory>
#include <algorithm>

struct Line {
    int a;
    int b;

    Line(int a, int b) : a(a), b(b) {}

    bool intersects(const Line other) const {
        int a1 = a, a2 = other.a;
        int b1 = b, b2 = other.b;
        if (a1 > b1) std::swap(a1, b1);
        if (a2 > b2) std::swap(a2, b2);

        int overlap_start = std::max(a1, a2);
        int overlap_end = std::min(b1, b2);

        if (overlap_start < overlap_end) {
            return true;
        } else {
            return false;
        }
    }

    Line intersection(Line other) const {
        int a1 = a, a2 = other.a;
        int b1 = b, b2 = other.b;
        if (a1 > b1) std::swap(a1, b1);
        if (a2 > b2) std::swap(a2, b2);

        int overlap_start = std::max(a1, a2);
        int overlap_end = std::min(b1, b2);

        if (overlap_start < overlap_end) {
            return {overlap_start, overlap_end};
        } else {
            return {INT_MIN, INT_MIN};
        }
    }
    bool isEmpty() const {
        return a == INT_MIN && b == INT_MIN;
    }
    bool operator==(const Line& other) const {
        return a == other.a && b == other.b;
    }
};

class SegmentTree {
public:
    bool xFirst;
    SegmentTree(const Line& boundary, bool xFirst) : xFirst(xFirst) {
        root = std::make_unique<Node>(boundary);
    }

    bool insert(Line x, Line y, Edge e) {
        return root->insert(x, y, e, xFirst);
    }

    int count(const Line& r) {
        return root->count(r);
    }

    void retrieve(const Point& p, std::vector<Edge> &res) {
        root->retrieve(p, res, xFirst);
    }

    void construct() {
        root->construct();
    }
    
    void serialize(const std::string& filename) const {
        std::ofstream ofs(filename, std::ios::binary);
        if (ofs.is_open()) {
            ofs.write(reinterpret_cast<const char*>(&xFirst), sizeof(xFirst));
            root->serialize(ofs);
            ofs.close();
        }
    }

    void deserialize(const std::string& filename) {
        std::ifstream ifs(filename, std::ios::binary);
        if (ifs.is_open()) {
            ifs.read(reinterpret_cast<char*>(&xFirst), sizeof(xFirst));
            root = std::make_unique<Node>(Line(0, 0)); // root node
            root->deserialize(ifs);
            ifs.close();
        }
    }

    void clear() {
        root->clear();
    }

private:
    struct Node {
        // t1, t2
        Line boundary;
        std::vector<std::pair<Line, Edge>> bound_edges;
        // s, t
        std::vector<bint> op_y1_;
        std::vector<std::vector<bint>> op_y2_;
        std::vector<std::vector<std::vector<Edge>>> op_values_;
        std::unique_ptr<Node> leftChild;
        std::unique_ptr<Node> rightChild;
        bool childEmpty = true;

        Node(const Line& boundary)
                : boundary(boundary) {}

        bool insert(Line x, Line y, Edge e, bool xFirst) {
            Line t_line(0, 0);
            if (xFirst) {
                t_line = x;
            } else {
                t_line = y;
            }
            auto sub_line = t_line.intersection(boundary);
            // non-overlap
            if (sub_line.isEmpty()) {
                return false;
            }
            // node is covered by the line
            if (sub_line == boundary) {
                // insert another line into this node
                if (xFirst) {
                    bound_edges.emplace_back(y, e);
                } else {
                    bound_edges.emplace_back(x, e);
                }
            } else { // some child nodes are covered by the rectangle
                if (childEmpty) {
                    split();
                }
                auto quadIndex = quadrantOf(t_line);
                for (auto &childIndex : quadIndex) {
                    if (childIndex == 0) {
                        leftChild->insert(x, y, e, xFirst);
                    } else {
                        rightChild->insert(x, y, e, xFirst);
                    }
                }
            }

            return true;
        }

        int count(Line r) {
            auto sub_r = r.intersection(boundary);
            // non-overlap
            if (sub_r.isEmpty()) {
                return 0;
            }
            // node is covered by the rectangle
            if (sub_r == boundary) {
                return 1;
            } else { // some child nodes are covered by the rectangle
                if (childEmpty) {
                    split();
                }
                int sum = 0;
                auto quadIndex = quadrantOf(r);
                for (auto &childIndex : quadIndex) {
                    if (childIndex == 0) {
                        sum += leftChild->count(r);
                    } else {
                        sum += rightChild->count(r);
                    }
                }
                return sum;
            }
        }

        void split() {
            int mid = (boundary.a + boundary.b) / 2;

            leftChild = std::make_unique<Node>(Line(boundary.a, mid));
            rightChild = std::make_unique<Node>(Line(mid, boundary.b));
            childEmpty = false;
        }

        void construct() {
            std::sort(bound_edges.begin(), bound_edges.end(), [](std::pair<Line, Edge> a, std::pair<Line, Edge> b) {
                if (a.first.a != b.first.a) {
                    return a.first.a < b.first.a;
                }

                return a.first.b > b.first.b;
            });

            if (!bound_edges.empty()) {
                bint target_y1 = bound_edges[0].first.a;
                bint target_y2 = bound_edges[0].first.b;
                Edge e = bound_edges[0].second;
                size_t y1_size = 0;
                size_t y2_size;
                op_y1_.push_back(target_y1);
                y1_size++;
                op_y2_.resize(y1_size);
                op_y2_[y1_size-1].push_back(target_y2);
                y2_size = 1;
                op_values_.resize(y1_size);
                op_values_[y1_size-1].resize(y2_size);
                op_values_[y1_size-1][y2_size-1].push_back(e);

                for (auto i = 1; i < bound_edges.size(); i++) {
                    bint y1 = bound_edges[i].first.a;
                    bint y2 = bound_edges[i].first.b;
                    e = bound_edges[i].second;

                    if (y1 != target_y1) {
                        op_y1_.push_back(y1);
                        y1_size++;
                        op_y2_.resize(y1_size);
                        op_y2_[y1_size-1].push_back(y2);
                        y2_size = 1;
                        op_values_.resize(y1_size);
                        op_values_[y1_size-1].resize(y2_size);
                        op_values_[y1_size-1][y2_size-1].push_back(e);
                        target_y1 = y1;
                        target_y2 = y2;
                        continue;
                    }

                    if (y2 != target_y2) {
                        op_y2_[y1_size-1].push_back(y2);
                        y2_size++;
                        op_values_[y1_size-1].resize(y2_size);
                        op_values_[y1_size-1][y2_size-1].push_back(e);
                        target_y2 = y2;
                        continue;
                    }

                    op_values_[y1_size-1][y2_size-1].push_back(e);
                }

                bound_edges.clear();
                bound_edges.shrink_to_fit();
            }

            if (!childEmpty) {
                if (!leftChild->bound_edges.empty()) leftChild->construct();
                if (!rightChild->bound_edges.empty()) rightChild->construct();
            }
        }

        [[nodiscard]] std::vector<int> quadrantOf(const Line& line) const {
            std::vector<int> res;

            if (line.intersects(leftChild->boundary)) {
                res.push_back(0);
            } if (line.intersects(rightChild->boundary)) {
                res.push_back(1);
            }

            return res;
        }

        void retrieve(const Point& p, std::vector<Edge> &resultSet, bool xFirst) const {
            for (int i = 0; i < op_y1_.size(); i++) {
                auto y1 = op_y1_[i];
                if (xFirst) {
                    if (p.y <= y1) break;
                } else {
                    if (p.y < y1) break;
                }
                for (int j = 0; j < op_y2_[i].size(); j++) {
                    auto y2 = op_y2_[i][j];
                    if (xFirst) {
                        if (y2 < p.y) break;
                    } else {
                        if (y2 <= p.y) break;
                    }
                    resultSet.insert(resultSet.end(), op_values_[i][j].begin(), op_values_[i][j].end());
                }
            }

            if (!childEmpty) {
                if (xFirst) {
                    if (leftChild->boundary.a <= p.x && p.x < leftChild->boundary.b) {
                        leftChild->retrieve(p, resultSet, xFirst);
                    } else {
                        rightChild->retrieve(p, resultSet, xFirst);
                    }
                } else {
                    if (leftChild->boundary.a < p.x && p.x <= leftChild->boundary.b) {
                        leftChild->retrieve(p, resultSet, xFirst);
                    } else {
                        rightChild->retrieve(p, resultSet, xFirst);
                    }
                }
            }
        }

        void serialize(std::ofstream& ofs) const {
            ofs.write(reinterpret_cast<const char*>(&boundary), sizeof(boundary));

            int y1_size = op_y1_.size();
            ofs.write(reinterpret_cast<const char*>(&y1_size), sizeof(int));
            ofs.write(reinterpret_cast<const char*>(op_y1_.data()), y1_size * sizeof(bint));
            for (int i = 0; i < y1_size; i++) {
                int y2_size = op_y2_[i].size();
                ofs.write(reinterpret_cast<const char*>(&y2_size), sizeof(int));
                ofs.write(reinterpret_cast<const char*>(op_y2_[i].data()), y2_size * sizeof(bint));
                for (int j = 0; j < y2_size; j++) {
                    int e_size = op_values_[i][j].size();
                    ofs.write(reinterpret_cast<const char*>(&e_size), sizeof(int));
                    ofs.write(reinterpret_cast<const char*>(op_values_[i][j].data()), e_size * sizeof(Edge));
                }
            }

            ofs.write(reinterpret_cast<const char*>(&childEmpty), sizeof(bool));
            if (!childEmpty) {
                leftChild->serialize(ofs);
                rightChild->serialize(ofs);
            }
        }

        void deserialize(std::ifstream& ifs) {
            ifs.read(reinterpret_cast<char*>(&boundary), sizeof(boundary));

            int y1_size;
            ifs.read(reinterpret_cast<char*>(&y1_size), sizeof(int));
            op_y1_.resize(y1_size);
            ifs.read(reinterpret_cast<char*>(op_y1_.data()), y1_size * sizeof(bint));
            op_y2_.resize(y1_size);
            op_values_.resize(y1_size);
            for (int i = 0; i < y1_size; i++) {
                int y2_size;
                ifs.read(reinterpret_cast<char*>(&y2_size), sizeof(int));
                op_y2_[i].resize(y2_size);
                op_values_[i].resize(y2_size);
                ifs.read(reinterpret_cast<char*>(op_y2_[i].data()), y2_size * sizeof(bint));
                for (int j = 0; j < y2_size; j++) {
                    int e_size;
                    ifs.read(reinterpret_cast<char*>(&e_size), sizeof(int));
                    op_values_[i][j].resize(e_size);
                    ifs.read(reinterpret_cast<char*>(op_values_[i][j].data()), e_size * sizeof(Edge));
                }
            }

            ifs.read(reinterpret_cast<char*>(&childEmpty), sizeof(bool));
            if (!childEmpty) {
                leftChild = std::make_unique<Node>(Line(0, 0));
                leftChild->deserialize(ifs);
                rightChild = std::make_unique<Node>(Line(0, 0));
                rightChild->deserialize(ifs);
            }
        }

        void clear() {
            op_y1_.clear();
            op_y1_.shrink_to_fit();
            op_y2_.clear();
            op_y2_.shrink_to_fit();
            op_values_.clear();
            op_values_.shrink_to_fit();
            if (!childEmpty) {
                leftChild->clear();
                rightChild->clear();
            }
        }
    };

    std::unique_ptr<Node> root;
};

#endif //HYBRID_H
