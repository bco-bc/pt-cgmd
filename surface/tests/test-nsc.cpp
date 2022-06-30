/*
 * Author: André H. Juffer.
 * Created on 18/05/2022, 19:03.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/surface/dotted-surface-generator.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include <iostream>
#include <fstream>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(simploce::util::Logger::LOGTRACE);

    std::vector<position_t> r;
    std::vector<radius_t> radii;

    radii.emplace_back(radius_t{0.15});
    radii.emplace_back(radius_t{0.15});

    r.emplace_back(position_t{0.0, 0.0, 0.0});
    r.emplace_back(position_t{0.2, 0.0, 0.0});

    auto result = simploce::dotted_surface_generator::general(r, radii);
    auto dots = result.first;
    std::clog << "Number of dots: " << dots.size() << std::endl;
    std::ofstream ostream;
    util::open_output_file(ostream, "/wrk3/tests/dotted-surface.dat");
    ostream << dots.size() << std::endl;
    for (auto& dot: dots ) {
        ostream << dot << std::endl;
    }
    ostream.close();

    return 0;
}