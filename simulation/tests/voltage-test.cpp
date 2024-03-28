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
#include "simploce/util/file.hpp"
#include "simploce/util/logger.hpp"
#include <cstdlib>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGDEBUG);

    std::ofstream ostream;
    std::string fn = "/wrk3/tests/voltage.dat";
    util::open_output_file(ostream, fn);

    //auto catalog = factory::particleSpecCatalog("/wrk3/tests/particle-specs.dat");
    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto box = factory::box(20.0);
    auto bc = factory::pbc1dSR(box, Direction::Z);
    bool mesoscale{false};
    dist_t dW = 3.0e-06;                          // In m.
    dist_t distance{20.0};                        // Box length in the z-direction.
    Voltage potential{el_field_t{0.0, 0.0, 1000.0}, bc, 78.5, mesoscale};
    auto particle = Atom::create("Na+", 0, "12345-98765", catalog->lookup("Na+"));
    std::cout << "Particle reset: " << particle->charge() << std::endl;
    real_t dr = distance() / 100.0;
    int n = int(distance() / dr);
    std:: cout << "n = " << n << std::endl;
    std::cout << "dr = " << dr << std::endl;
    ostream.setf(std::ios::scientific);
    position_t r{0.0, 0.0, 0.0};
    for (int i = -(n-1); i != n; ++i) {
        r[2] = i * dr;
        particle->position(r);
        auto result = potential.operator()(particle);  // Energy and force on particle.
        ostream << r[2] << ' ' << result.first << " " << result.second[2] << std::endl;
    }
    ostream.close();
    std::cout << "Data written to '" << fn << "'." << std::endl;

    return EXIT_SUCCESS;
}

