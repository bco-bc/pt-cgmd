//
// Created by ajuffer on 6/2/22.
//

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/util/logger.hpp"
#include <iostream>
#include <string>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGINFO);

    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto factory = factory::particleSystemFactory(catalog);

    std::string fn = "/wrk3/tests/1CUS.pdb";
    auto pdb = factory->fromPDB(fn);

    std::cout << "Number of atoms: " << pdb->numberOfParticles() << std::endl;
    std::cout << "Number of atom groups: " << pdb->numberOfParticleGroups() << std::endl;
    std::cout << "Number of free atoms: " << pdb->numberOfFreeParticles() << std::endl;

    std::clog << *pdb << std::endl;

    return 0;
}