/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */


#ifndef SURFACE_DOTTED_SURFACE_GENERATOR_HPP
#define SURFACE_DOTTED_SURFACE_GENERATOR_HPP

#include "./types/srf-types.hpp"
#include <vector>
#include <utility>

namespace simploce {

    struct dotted_surface_generator {

        /**
         * Returns dotted surface of sphere.
         * @param radius Sphere radius.
         * @param density Number of dots per area.
         * @return Dots.
         */
        static std::pair<std::vector<dot_t>, area_t> spherical(radius_t radius, real_t density);

        /**
         * Returns a dotted surface of collection of spherical objects.
         * @param positions Objects' positions.
         * @param radii Objects' radii. Must be non-zero values.
         * @return Dots.
         */
        static std::pair<std::vector<dot_t>, area_t> general(std::vector<position_t>& positions,
                                                             std::vector<radius_t>& radii);
    };
}

#endif //SURFACE_DOTTED_SURFACE_GENERATOR_HPP
