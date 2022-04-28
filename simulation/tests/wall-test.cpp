/*
 * Author: André H. Juffer.
 * Created on 22/12/2021, 16:13.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/wall.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include <cstdlib>

using namespace simploce;

int main() {
    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto box = factory::box(6.0);
    auto bc = factory::oneDimensionBoundaryCondition(box, Direction::Z);
    real_t C12 = 3.0e-08;
    real_t C6 = 0.00005;
    dist_t distanceToPlane{7.0};
    int coordinate = 1;
    FlatSurface flatSurface1{Plane::ZX, distanceToPlane};
    //std::cout << "Flat surface1: " << flatSurface1.toString() << std::endl;
    Wall wall1(C12, C6, bc, flatSurface1);
    FlatSurface flatSurface2{Plane::ZX, dist_t {0.0}};
    //std::cout << "Flat surface2: " << flatSurface2.toString() << std::endl;
    Wall wall2(C12, C6, bc, flatSurface2);

    auto particle = Atom::create("1345x", 0, "test", catalog->lookup("Na+"));
    real_t dr = 0.01;
    int n = distanceToPlane() / dr;
    std::cout.setf(std::ios::scientific);

    position_t r{1.0, dr, 1.0};
    particle->position(r);
    //auto result = wall1.operator()(particle);
    // std::cout << r << ' ' << result.first << result.second << std::endl;

    for (int i = 1; i != n; ++i) {
        r[coordinate] = i * dr;
        particle->position(r);
        auto result1 = wall1.operator()(particle);
        auto result2 = wall2.operator()(particle);
        auto energy = result1.first + result2.first;
        auto f = result1.second + result2.second;
        //auto energy = result2.first;
        //auto f = result2.second;
        std::cout << r << ' ' << energy << ' ' << f << std::endl;
    }

    return EXIT_SUCCESS;
}