/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#ifndef SURFACE_TRIANGULATION_HPP
#define SURFACE_TRIANGULATION_HPP

#include "polyhedron-generator.hpp"
#include <iostream>

namespace simploce {

    struct triangulation : public polyhedron_generator {

        polyhedron_ptr_t cubic(length_t sideLength) override;

        polyhedron_ptr_t spherical(radius_t radius, std::size_t numberOfTriangles) override;

        polyhedron_ptr_t parse(std::istream& stream) override;

    };
}
#endif //SURFACE_TRIANGULATION_HPP
