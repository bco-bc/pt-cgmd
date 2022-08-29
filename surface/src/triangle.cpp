/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#include "simploce/surface/triangle.hpp"
#include "simploce/surface/edge.hpp"
#include "simploce/surface/vertex.hpp"
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simploce {

    face_ptr_t
    Triangle::create(const vertex_ptr_t &v1,
                     const vertex_ptr_t &v2,
                     const vertex_ptr_t &v3) {
        std::vector<vertex_ptr_t> vertices;
        vertices.emplace_back(v1);
        vertices.emplace_back(v2);
        vertices.emplace_back(v3);
        std::vector<edge_ptr_t> edges;
        edge_ptr_t edge = Edge::create(v1, v2);
        edges.emplace_back(edge);
        edge = Edge::create(v2, v3);
        edges.emplace_back(edge);
        edge = Edge::create(v3, v1);
        edges.emplace_back(edge);
        auto face = std::make_shared<Triangle>(vertices, edges);
        face->validate_();
        return face;
    }

    Triangle::Triangle(std::vector<vertex_ptr_t> vertices, std::vector<edge_ptr_t> edges) :
        Face(std::move(vertices), std::move(edges)) {
    }

    area_t Triangle::area() {
        if (!hasArea_) {
            auto edges = this->edges();
            auto a = edges[0]->length();
            auto b = edges[1]->length();
            auto c = edges[2]->length();
            auto s = 0.5 * (a + b + c);
            auto t = s * (s - a) * (s - b) * (s - c);
            area_ = std::sqrt(t);
            hasArea_ = true;
        }
        return area_;
    }

    void Triangle::validate_() {
        auto vertices = this->vertices();
        auto edges = this->edges();
        if (vertices.size() != 3) {
            throw std::domain_error("Triangle: Must have exactly 3 vertices.");
        }
        if (edges.size() != 3) {
            throw std::domain_error("Triangle: Must have exactly 3 edges.");
        }

        if ( vertices[0] == vertices[1] ||
             vertices[1] == vertices[2] ||
             vertices[2] == vertices[0] ) {
            throw std::domain_error("Triangle: 2 or 3 vertices are identical.");
        }
        int n = 0;

        for (auto& edge: edges) {
            auto start = edge->start();
            for (auto& vertex: vertices) {
                if (start == vertex) {
                    n++;
                }
            }
        }
        if (n != 3)
            throw std::domain_error("A triangle's edges must refer to its own vertices");
    }

 }
