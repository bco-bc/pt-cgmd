/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/logger.hpp"
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fnSpecs = "/localdisk/resources/particles-specs.dat";
    // std::string fnSpecs = "/wrk3/simulation/NaClNextToChargedSurface/no-external-field/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fnSpecs);
    auto catalog = factory::particleSpecCatalog(stream);
    stream.close();

    std::string fnInteractions = "/localdisk/resources/NaCl-uniform-srf-cg-density-virtual-planes.dat";
    //std::string fnInteractions = "/localdisk/resources/interaction-parameters.dat";
    //std::string fnInteractions = "/wrk3/simulation/NaClNextToChargedSurface/no-external-field/NaCl-interaction-parameters.dat";
    util::open_input_file(stream, fnInteractions);
    auto forceField = factory::forceField(stream, catalog);
    stream.close();

    std::cout << "Force field:" << std::endl;
    std::cout << *forceField << std::endl;
    std::cout << std::endl;

    auto spec1 = catalog->lookup("Na+");
    auto spec2 = catalog->lookup("Cl-");
    auto param = forceField->hardSphereShiftedForce(spec1, spec2);
    //auto param = forceField->hardSphereLekner(spec1, spec2);
    std::cout << "Na+ - Cl- interaction param: " << param << std::endl;

    // Should cause exception.
    forceField->uniformSurfaceChargeDensity();

    return EXIT_SUCCESS;
}