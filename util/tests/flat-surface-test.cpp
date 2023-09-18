/*
 * Author: André H. Juffer.
 * Created on 21/12/2021, 15:51.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/flat-surface.hpp"
#include <iostream>

using namespace simploce;

int main() {

    FlatSurface flatSurface1{Plane::ZX, dist_t {0.0}};
    std::cout << "Surface #1: " << flatSurface1.toString() << std::endl;
    FlatSurface flatSurface2{Plane::ZX, dist_t {7.0}};
    std::cout << "Surface #2: " << flatSurface2.toString() << std::endl;

    position_t r1{1, 1, 1.1};
    position_t r2{0.4, 0.89, -1.1};

    auto pair = flatSurface1.distanceTo(r1);
    std::cout << "Distance r1 (" << r1 << " ) to surface 1: " << pair.first << ",  distance vector: " << pair.second << std::endl;
    std::cout << "Answer: distance = 1.0, distance vector: 0.0, 1.0, 0.0" << std::endl;

    pair = flatSurface2.distanceTo(r1);
    std::cout << "Distance r1 (" << r1 << " ) to surface 2: " << pair.first << ",  distance vector: " << pair.second << std::endl;
    std::cout << "Answer: distance = 6.0, distance vector: 0.0, -6.0, 0.0" << std::endl;

    auto normal = flatSurface2.unitVectorPerpendicularTo();
    std::cout << "Normal vector: " << normal << std::endl;


    return EXIT_SUCCESS;
}