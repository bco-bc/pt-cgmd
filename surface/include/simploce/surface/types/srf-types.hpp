/*
 * Types for surface calculations.
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#ifndef SURFACE_SRF_TYPES_HPP
#define SURFACE_SRF_TYPES_HPP

#include "simploce/types/u-types.hpp"
#include <memory>

namespace simploce {

    // Forward declarations.
    class Vertex;
    class Edge;
    class Face;
    class Polyhedron;

    // Dot: Point of a dotted surface.
    using dot_t = position_t;

    using vertex_ptr_t = std::shared_ptr<Vertex>;

    using edge_ptr_t = std::shared_ptr<Edge>;

    using face_ptr_t = std::shared_ptr<Face>;

    using polyhedron_ptr_t = std::shared_ptr<Polyhedron>;

}

#endif //SURFACE_SRF_TYPES_HPP
