//
// Created by ajuffer on 5/19/22.
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
    auto box = factory::box(10, 10, 40);
    //auto box = factory::box(10, 10, 10);

    auto polymerSolution =
            factory->polymerSolution(box,
                                     15,
                                     "PMU",
                                     400,
                                     1.0,
                                     4000,
                                     "H2Om",
                                     1.0,
                                     true);

/*
    auto polymerSolution =
            factory->polymerSolution(box,
                                     15,
                                     "PMU",
                                     0,
                                     1.0,
                                     12000,
                                     "H2Om",
                                     2.0,
                                     true);
*/
    std::ofstream stream;
    util::open_output_file(stream, "/wrk3/tests/polymer-solution.ps");
    //util::open_output_file(stream, "/wrk3/tests/water-dpd.ps");
    stream << *polymerSolution << std::endl;
    stream.close();
}