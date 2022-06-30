/*
 * Author: André H. Juffer.
 * Created on 16/06/2022.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/flat-triangles-calculator.hpp"
#include "simploce/surface/triangulation.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/logger.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    triangulation triangulator;
    auto surface = triangulator.spherical(2.0, 240);
    auto param = std::make_shared<param_t>();
    param->put<real_t>("bem.solvent.ka", 0.06);
    param->put<real_t>("bem.solvent.eps", 80.0);
    param->put<real_t>("bem.solute.eps", 20.0);

    bem_calc_ptr_t calculator = FlatTrianglesCalculator::create(param, surface);

    calculator->surfaceMatrix();

    std::vector<position_t> positions;
    positions.emplace_back(position_t{});
    std::vector<charge_t> charges;
    charges.emplace_back(charge_t{1.0});
    calculator->rightHandSide(positions, charges);

    calculator->solve();

    std::vector<position_t> points;
    points.emplace_back(position_t{});
    auto potentials = calculator->electricPotentials(points);

    for (std::size_t i = 0; i != points.size(); ++i) {
        std::cout << points[i] << " " << potentials[i] << std::endl;
    }

    return EXIT_SUCCESS;

}