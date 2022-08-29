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

    std::pair<position_t, normal_t>
    Face::center() {
        if (!hasCenter_) {
            auto normal = this->normal();
            position_t midPoint{0.0, 0.0, 0.0};
            for (auto &vertex: vertices_) {
                midPoint += vertex->position();
            }
            midPoint /= vertices_.size();
            center_ = std::make_pair(midPoint, normal);
            hasCenter_ = true;
        }
        return center_;
    }

    normal_t
    Face::normal() {
        if ( !hasNormal_) {
            normal_ = {0.0, 0.0, 0.0};
            for (auto &vertex: vertices_) {
                normal_ += vertex->normal();
            }
            normal_ /= vertices_.size();
            normal_ /= norm<real_t>(normal_);
            normal_ = normal_t{normal_};
            hasNormal_ = true;
        }
        return normal_;
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
