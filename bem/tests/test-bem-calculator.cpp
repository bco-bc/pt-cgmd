/*
 * Author: André H. Juffer.
 * Created on 16/06/2022.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/flat-tri-nodes-vertices-calculator.hpp"
#include "simploce/bem/flat-tri-nodes-tri-calculator.hpp"
#include "simploce/bem/curved-tri-nodes-vertices-calculator.hpp"
#include "simploce/surface/triangulation.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    radius_t radius{3.0};
    triangulation triangulator;
    auto surface = triangulator.spherical(radius, 60);
    std::string fn = "/wrk3/tests/triangulated-surface.dat";
    std::ofstream stream;
    util::open_output_file(stream, fn);
    stream << *surface << std::endl;
    stream.close();

    auto param = std::make_shared<param_t>();
    //param->put<real_t>("bem.solvent.ka", 0.06);
    param->put<real_t>("bem.solvent.ka", 0.0);
    param->put<real_t>("bem.solvent.eps", 78.5);
    param->put<real_t>("bem.solute.eps", 1.0);

    //auto calculator = FlatTriNodesVerticesCalculator::create(param, surface);
    //auto calculator = FlatTriNodesTriCalculator::create(param, surface);
    auto calculator = CurvedTriNodesVerticesCalculator::create(param, surface);

    calculator->surfaceMatrix();

    fn = "/wrk3/tests/reaction-potentials-bem.dat";
    util::open_output_file(stream, fn);
    std::vector<position_t> positions;
    positions.emplace_back(position_t{});
    std::vector<charge_t> charges;
    charges.emplace_back(charge_t{1.0});
    real_t dz = 0.01;
    int n = int(radius()/dz);

    for (int i = 0; i != n; ++i) {
        auto v = i * dz;
        positions[0] = {0, 0, v};
        auto d = norm<real_t>(positions[0]);
        if ( std::fabs(d - radius()) > 1.0e-05) {
            calculator->rightHandSide(positions, charges);
            calculator->solve();
            std::vector<el_pot_t> potentials;
            if (d < radius()) {
                potentials = calculator->reactionPotentialSolute(positions);
            } else {
                potentials = calculator->reactionPotentialSolvent(positions);
            }
            auto rp = potentials[0] / units::mu<real_t>::V_to_kJ_mol_e;

            std::cout << d << rp << std::endl;
            stream << d << rp << std::endl;
        }
    }
    stream.close();

    return EXIT_SUCCESS;

}