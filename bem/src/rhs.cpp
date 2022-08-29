/*
 * Author: André H. Juffer.
 * Created on 12/07/2022.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/rhs.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/units/units-mu.hpp"

namespace simploce {
    namespace rhs {

        BEMData::vector_t rightHandSide(std::vector<BEMData::node_t>& nodes,
                                        std::vector<position_t> &rcg,
                                        std::vector<charge_t> &cg,
                                        real_t epsSolvent,
                                        real_t epsSolute) {
            util::Logger logger("simploce::rhs::rightHandSide()");
            logger.trace("Entering.");

            if (nodes.empty()) {
                util::logAndThrow(logger, "Nodes not provided.");
            }
            if (rcg.size() != cg.size()) {
                util::logAndThrow(logger, "Length cg and positions vector are not the same.");
            }
            if (rcg.empty() || cg.empty()) {
                util::logAndThrow(logger, "Positions and/or cg not provided.");
            }

            auto epsRatio = epsSolvent / epsSolute;
            BEMData::vector_t b{nodes.size()};
            b.setZero();
            auto FOUR_PI_E0 = units::mu<real_t>::FOUR_PI_E0;
            auto f1 = 2.0 / ((1.0 + epsRatio) * FOUR_PI_E0 * epsSolute);
            for (auto i = 0; i != cg.size(); ++i) {
                auto Q = cg[i]();
                auto rQ = rcg[i];
                for (auto node: nodes) {
                    auto& r0 = node.r;
                    auto& index = node.index;
                    auto R = rQ - r0;
                    auto dis = norm<real_t>(R);
                    b[index] += f1 * Q / dis;
                }
            }

            //std::clog << "RHS: "<< std::endl << b << std::endl;

            logger.trace("Leaving.");
            return std::move(b);
        }
    }
}