/*
 * Author: André H. Juffer.
 * Created on 30/11/2021, 16:09.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/voltage.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>

using namespace simploce;

int main() {
    std::ofstream ostream;
    std::string fn = "/wrk3/tests/voltage.dat";
    util::open_output_file(ostream, fn);

    std::cout << "V to kJ/(mol e): " << units::mu<real_t>::V_to_kJ_mol_e << std::endl;
    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto box = factory::box(6.0);  // nm.
    auto bc = factory::pbc(box);
    auto direction = Direction::valueOf('z');
    voltage_t voltage(0.050); // in V, or 50 mV
    dist_t distance{6.0};
    real_t eps_r = 1.0;
    Voltage potential{voltage, distance, eps_r, bc, direction};

    auto particle = Atom::create("Na+", 0, "12345-98765", catalog->lookup("Na+"));
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
        auto result = potential.operator()(particle);  // Electric potential and field.
        std::cout << r << ' ' << result.first << result.second << std::endl;
        ostream << r << ' ' << result.first << result.second << std::endl;
    }
    ostream.close();
    std::cout << "Data written to '" << fn << "'." << std::endl;

    return EXIT_SUCCESS;
}

