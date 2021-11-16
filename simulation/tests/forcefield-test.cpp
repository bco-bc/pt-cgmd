/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on on 11/8/21.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/force-field.hpp"
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
    std::ifstream stream;
    util::open_input_file(stream, fnSpecs);
    auto catalog = factory::particleSpecCatalog(stream);
    stream.close();

    std::string fnInteractions = "/localdisk/resources/interaction-parameters.dat";
    util::open_input_file(stream, fnInteractions);
    auto forceField = factory::obtainFrom(stream, catalog);
    stream.close();

    std::cout << "Force field:" << std::endl;
    std::cout << *forceField << std::endl;
    std::cout << std::endl;

    forceField->relativePermittivityInsideCutoff(4.0);
    std::cout << "Relative permittivity: " << forceField->relativePermittivityInsideCutoff() << std::endl;

    return EXIT_SUCCESS;
}