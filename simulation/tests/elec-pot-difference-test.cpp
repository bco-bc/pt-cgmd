/*
 * Author: André H. Juffer.
 * Created on 30/11/2021, 16:09.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/elec-pot-difference.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/units/units-mu.hpp"
#include <cstdlib>

using namespace simploce;

int main() {
    std::cout << "V to kJ/(mol e): " << units::mu<real_t>::V_to_kJ_mol_e << std::endl;
    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto box = factory::box(6.0);
    auto bc = factory::boundaryCondition(box);
    auto direction = Direction::valueOf('y');
    ElectricPotentialDifference epd{0.050, 6.3, 78.5, bc, direction};

    auto particle = Atom::create("1345x", 0, "test", catalog->lookup("Na+"));
    real_t dr = 0.01;
    int n = 41;
    std::cout.setf(std::ios::scientific);
    position_t r{0.0, 0.0, 0.0};
    for (int i = 0; i != n; ++i) {
        if (direction == Direction::Z) {
            r[2] = i * dr;
        } else if (direction == Direction::X) {
            r[0] = i * dr;
        } else {
            r[1] = i * dr;
        }
        particle->position(r);
        auto result = epd.operator()(particle);
        std::cout << r << ' ' << result.first << result.second << std::endl;
    }

    return EXIT_SUCCESS;
}

