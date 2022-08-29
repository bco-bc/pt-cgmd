/*
 * Author: André H. Juffer.
 * Created on 06/06/22, 21:39.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/surface/pentagon.hpp"
#include "simploce/surface/edge.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/math-constants.hpp"
#include <memory>
#include <utility>

namespace simploce {

    face_ptr_t
    Pentagon::create(const vertex_ptr_t &v1,
                     const vertex_ptr_t &v2,
                     const vertex_ptr_t &v3,
                     const vertex_ptr_t &v4,
                     const vertex_ptr_t &v5) {
        std::vector<vertex_ptr_t> vertices;
        vertices.emplace_back(v1);
        vertices.emplace_back(v2);
        vertices.emplace_back(v3);
        vertices.emplace_back(v4);
        vertices.emplace_back(v5);
        std::vector<edge_ptr_t> edges;
        edge_ptr_t edge = Edge::create(v1, v2);
        edges.emplace_back(edge);
        edge = Edge::create(v2, v3);
        edges.emplace_back(edge);
        edge = Edge::create(v3, v4);
        edges.emplace_back(edge);
        edge = Edge::create(v4, v5);
        edges.emplace_back(edge);
        edge = Edge::create(v5, v1);
        edges.emplace_back(edge);
        auto face = std::make_shared<Pentagon>(vertices, edges);
        face->validate_();
        return face;

    }

    Pentagon::Pentagon(std::vector<vertex_ptr_t> vertices, std::vector<edge_ptr_t> edges) :
            Face(std::move(vertices), std::move(edges)) {
    }

    area_t
    Pentagon::area() {
        static auto PI = math::constants<real_t>::PI;
        static util::Logger logger("simploce::surface::Pentagon::area()");
        logger.trace("Entering.");

        // Find average edge length.
        length_t total{0.0};
        auto edges = this->edges();
        for (auto& edge: edges) {
            auto length = edge->length();
            total += length;
        }
        auto average = total / real_t(edges.size());
        logger.debug(util::toString(total) + ": Average side (edge) length.");
        auto t1 = std::tan(3.0 * PI / 10.0);
        area_t area = 5.0 * average() * average() * std::tan(3.0 * PI / 10.0) / 4.0;
        logger.trace("Leaving.");
        return area;
    }

    void
    Pentagon::validate_() {
        auto vertices = this->vertices();
        auto edges = this->edges();
        if (vertices.size() != 5) {
            throw std::domain_error("Pentagon: Must have exactly 5 vertices.");
        }
        if (edges.size() != 5) {
            throw std::domain_error("Pentagon: Must have exactly 5 edges.");
        }
        bool equal{false};
        for (std::size_t i = 0; i != 4; ++i) {
            for (std::size_t j = i + 1; j != 5; ++j) {
                equal = vertices[i] == vertices[j];
            }
        }
        if (equal) {
            throw std::domain_error("Pentagon: 2 or more vertices are identical.");
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
        if (n != 5)
            throw std::domain_error("A pentagon's edges must refer to its own vertices");
    }
}
