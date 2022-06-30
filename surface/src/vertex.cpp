/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/17/22.
 */

#include "simploce/surface/vertex.hpp"
#include <memory>

namespace simploce {

    vertex_ptr_t Vertex::create(const position_t& r, const normal_t& un) {
        return std::make_shared<Vertex>(r, un);
    }

    Vertex::Vertex(const position_t r, const normal_t un) :
        index_{INDEX_++}, r_{r}, un_{un} {
    }

    int Vertex::index() const {
        return index_;
    }

    real_t Vertex::operator [] (int k) const {
        return r_[k];
    }

    position_t Vertex::position() const {
        return r_;
    }

    normal_t Vertex::normal() const {
        return un_;
    }

    int Vertex::INDEX_ = 0;

    void Vertex::position_(const position_t& r) {
        r_ = r;
    }

    void Vertex::normal_(const normal_t &normal) {
        un_ = normal;
    }
}