/*
 * Author: André H. Juffer.
 * Created on 19/12/2021, 13:35.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/sf.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <string>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    // Pair potential.
    std::string fileName = "/localdisk/resources/particles-specs.dat";
    auto catalog = factory::particleSpecCatalog(fileName);
    fileName = "/localdisk/resources/interaction-parameters.dat";
    auto forceField = factory::forceField(fileName, catalog);
    dist_t cutoff = 2.6;
    auto box = factory::box(2.0 * cutoff());
    auto bc = factory::boundaryCondition(box);
    SF pairPotential{cutoff, forceField, box, bc};

    // Particles.
    auto T1 = catalog->lookup("T1");
    auto T2 = catalog->lookup("T2");
    p_ptr_t p1 = Atom::create("1", 0, "p1", T1);
    p1->position(position_t{0.0, 0.0, 0.0});
    p_ptr_t p2 = Atom::create("2", 1, "p2", T2);
    p2->position(position_t{0.0, 0.0, 0.0});

    // Calculate potential.
    real_t dz = 0.1;
    int n = util::nint(cutoff() / dz);
    for (int i = 1; i < n; ++i) {
        real_t z = i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = pairPotential.operator()(p1, p2);
        std::cout << z << conf::SPACE << result.first << conf::SPACE << result.second << std::endl;
    }


}