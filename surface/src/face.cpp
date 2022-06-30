/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#include "simploce/surface/face.hpp"
#include "simploce/surface/edge.hpp"
#include "simploce/surface/vertex.hpp"

namespace simploce {

    Face::~Face() = default;

    Face::Face(std::vector<vertex_ptr_t> vertices, std::vector<edge_ptr_t> edges) :
        vertices_{std::move(vertices)}, edges_{std::move(edges)} {
        for (const auto& vertex: vertices_) {
            indices_.push_back(vertex->index());
        }
    }

    std::vector<int>
    Face::indices() const {
        return indices_;
    }

    std::vector<vertex_ptr_t>
    Face::vertices() const {
        return vertices_;
    }

    std::vector<edge_ptr_t>
    Face::edges() const {
        return edges_;
    }

    void
    Face::resetEdges(const std::map<id_t, edge_ptr_t> &edges) {
       std::vector<id_t> ids{};
       for (auto& edge: edges_) {
           ids.emplace_back(edge->id());
       }
       edges_.clear();
       for (auto& id: ids) {
           auto edge = edges.at(id);
           edges_.emplace_back(edge);
       }
       this->validate_();
    }

}
