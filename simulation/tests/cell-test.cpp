/*
 * Author: André H. Juffer.
 * Created on 31/08/22, 20:44.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/cell.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

int main() {

    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    auto catalog =
            factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto factory = simploce::factory::particleSystemFactory(catalog);
    auto particleSystem = factory->polarizableWater();

    Cell::location_t location(0,0,0);
    position_t r(0.5, 0.5, 0.5);
    Cell cell(location, r);

    particleSystem->doWithAll<void>([&cell] (const std::vector<p_ptr_t>& all) {
        for (const auto& p : all) {
            cell.assign(p);
        }
    });

    std::cout << "Number of particles: " << particleSystem->numberOfParticles() << std::endl;
    std::cout << "Cell particles size: " << cell.particles().size() << std::endl;

    return EXIT_SUCCESS;
}
