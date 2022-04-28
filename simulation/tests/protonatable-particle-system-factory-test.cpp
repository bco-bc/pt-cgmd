/*
 * Author: André H. Juffer.
 * Created on 29/10/2021, 14:59.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/types/s-types.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include <cstdlib>

using namespace simploce;

void protonatableWater(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating coarse grained protonatable polarizable water model:" << std::endl;
    auto factory = factory::protonatableParticleSystemFactory(catalog);
    box_ptr_t box = factory::box(length_t{7.27});
    auto cg = factory->protonatablePolarizableWater(box);
    std::cout << *cg << std::endl;
    std::cout << std::endl;
}

void formicAcid(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating coarse grained HCOOH solution:" << std::endl;
    auto factory = factory::protonatableParticleSystemFactory(catalog);
    box_ptr_t box = factory::box(length_t{7.27});
    auto cg = factory->formicAcid(box);
    std::cout << *cg << std::endl;
    std::cout << std::endl;
}

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fileName =
            "/localdisk/resources/particles-specs.dat";
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fileName);
    std::cout << *catalog << std::endl << std::endl;

    //protonatableWater(catalog);
    formicAcid(catalog);

    return (EXIT_SUCCESS);
}


