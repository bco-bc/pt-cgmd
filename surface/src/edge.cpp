/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#include "simploce/surface/edge.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/types/cvector_t.hpp"
#include <stdexcept>
#include <utility>

namespace simploce {

    /**
     * Returns identifier for an edge between given start and end vertex. Identifier is based the
     * vertices' indexes.
     * @param start Start vertex.
     * @param end End vertex.
     * @return Identifier.
     */
    static id_t
    generateId_(const vertex_ptr_t& start, const vertex_ptr_t& end) {
        auto indexStart = start->index();
        auto indexEnd = end->index();
        auto id = indexStart > indexEnd ?
                std::to_string(indexEnd) + "-" + std::to_string(indexStart) :
                std::to_string(indexStart) + "-" + std::to_string(indexEnd);
        return std::move(id_t{id});
    }

    edge_ptr_t Edge::create(const vertex_ptr_t &start, const vertex_ptr_t &end) {
        return std::make_shared<Edge>(start, end);
    }

    Edge::Edge(vertex_ptr_t start, vertex_ptr_t end) :
            id_{}, start_{std::move(start)}, end_{std::move(end)} {
        this->validate_();
        id_ = generateId_(start_, end_);
    }

    id_t
    Edge::id() const {
        return id_;
    }

    bool
    Edge::operator == (const Edge &edge) {
        return ((start_ == edge.start_ && end_ == edge.end_) ||
                (start_ == edge.end_ && end_ == edge.start_) );
    }

    vertex_ptr_t
    Edge::start() const {
        return start_;
    }

    vertex_ptr_t
    Edge::end() const {
        return end_;
    }

    length_t
    Edge::length() const {
        return norm<real_t>(start_->position() - end_->position());
    }

    void Edge::validate_() {
        if ( !start_ ) {
            throw std::domain_error("Edge: Start vertex must be provided.");
        }
        if ( !end_ ) {
            throw std::domain_error("Edge: End vertex must be provided.");
        }
        if (start_ == end_) {
            throw std::domain_error("Edge: Start and end vertex must not be identical.");
        }
        auto dv = start_->position() - end_->position();
        if (norm_square<real_t>(dv) <= std::numeric_limits<real_t>::epsilon()) {
            throw std::domain_error("Edge: Start and end vertex must not be identical.");
        }

    }
}