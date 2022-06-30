/*
 * Author: André H. Juffer.
 * Created on 06/06/22, 22:31.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SURFACE_SPHERE_HPP
#define SURFACE_SPHERE_HPP

#include "simploce/types/u-types.hpp"
#include "simploce/surface/types/srf-types.hpp"
#include <utility>
#include <vector>

namespace simploce {
    namespace sphere {

        /**
         * Creates a spherical dodecahedron as a collection of pentagons.
         * @param radius Radius of sphere.
         * @return Vertices and pentagons. Edges are locally defined in pentagons only.
         */
        std::pair<std::vector<vertex_ptr_t>, std::vector<face_ptr_t>> dodecahedron(const radius_t& radius = 1.0);

        /**
         * Constructs 60 triangles from 12 pentagons on the surface of a sphere.
         * @param polyhedron Surface with 12 pentagons.
         * @param radius Radius of sphere.
         * @return 32 vertices, 60 triangles.
         */
        std::pair<std::vector<vertex_ptr_t>, std::vector<face_ptr_t>> triangles60(const polyhedron_ptr_t& polyhedron,
                                                                                  const radius_t& radius = 1.0);

        /**
         * Further divide triangles on the surface of a sphere. The division continues until the number of generated
         * triangles is equal to or larger than the requested number of triangles.
         * @param polyhedron Triangulated surface.
         * @param numberOfTriangles Requested number of triangles.
         * @param radius Radius of sphere.
         * @return Vertices, triangles.
         */
        std::pair<std::vector<vertex_ptr_t>, std::vector<face_ptr_t>> divide(const polyhedron_ptr_t& polyhedron,
                                                                             std::size_t numberOfTriangles = 960,
                                                                             const radius_t& radius = 1.0);
    }
}

#endif //SURFACE_SPHERE_HPP
