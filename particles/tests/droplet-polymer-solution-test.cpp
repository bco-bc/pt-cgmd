//
// Created by ajuffer on 8/29/22.
//

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include <iostream>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(simploce::util::Logger::LOGDEBUG);

    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto factory = simploce::factory::particleSystemFactory(catalog);
    auto box = factory::box(20, 20, 40);

    auto particleSystem =
            factory->dropletPolymerSolution(box,
                                            15,
                                            "PMU",
                                            0,
                                            1.0,
                                            24000,
                                            "H2Om",
                                            24000,
                                            "DROm",
                                            1.0);

    std::ofstream stream;
    util::open_output_file(stream, "/wrk3/tests/droplets-polymer-solution.ps");
    //util::open_output_file(stream, "/wrk3/tests/water-dpd.ps");
    stream << *particleSystem << std::endl;
    stream.close();

    return EXIT_SUCCESS;
}
