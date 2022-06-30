/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#ifndef SURFACE_POLYHEDRON_GENERATOR_HPP
#define SURFACE_POLYHEDRON_GENERATOR_HPP

#include "./types/srf-types.hpp"
#include <iostream>

namespace simploce {

    /**
     * Generates polyhedron for several shapes.
     */
    struct polyhedron_generator {

        virtual ~polyhedron_generator();

        /**
         * Returns cube.
         * @param sideLength Side length cube.
         * @return Polyhedron.
         */
        virtual polyhedron_ptr_t cubic(length_t sideLength) = 0;

        /**
         * Returns sphere.
         * @param radius Radius sphere.
         * @param numberOfTriangles Requested number of triangles.
         * @return Polyhedron.
         */
        virtual polyhedron_ptr_t spherical(radius_t radius, std::size_t numberOfTriangles) = 0;

        /**
         * Parses polyhedron from input stream.
         * @param stream Input stream.
         * @return Polyhedron.
         */
        virtual polyhedron_ptr_t parse(std::istream& stream) = 0;

        /**
         * Maps given polyhedron onto the given dotted surface. That is, all its vertices will
         * be placed on the dotted surface. The dot density (number of points per unit of area)
         * should be relatively high (> 500)
         * @param dots Dots. This must represent a closed and simple surface. May be modified upon return.
         * @param polyhedron Polyhedron.
         */
        static void mapOnto(std::vector<dot_t>& dots,
                            polyhedron_ptr_t& polyhedron);
    };
}
#endif //SURFACE_POLYHEDRON_GENERATOR_HPP
