/**
 * Author: André H. Juffer.
 * Created on 14/11/2023
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/halve-attractive-hp.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <limits>

namespace simploce {

    HalveAttractiveHP::HalveAttractiveHP(simploce::ff_ptr_t forceField, simploce::bc_ptr_t bc) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)} {
        util::Logger logger("simploce::HalveAttractiveHP::HalveAttractiveHP()");
        if (!forceField_) {
            util::logAndThrow(logger, "Missing force field.");
        }
        if (!bc_) {
            util::logAndThrow(logger, "Missing boundary conditions.");
        }
    }

    std::pair<energy_t, force_t>
    HalveAttractiveHP::operator()(const simploce::p_ptr_t &p1, const simploce::p_ptr_t &p2) {
        static util::Logger logger("simploce::HalveAttractiveHP::operator()()");
        logger.trace("Entering.");

        static auto eps = real_t(std::numeric_limits<float>::min());

        // Get r0 and fc.
        auto params = forceField_->halveAttractiveHarmonic(p1->spec(), p2->spec());
        auto r0 = params.first;
        auto fc = params.second;

        // Current positions.
        const auto &r1 = p1->position();
        const auto &r2 = p2->position();

        // Distance.
        auto rij = bc_->apply(r1, r2);
        auto Rij = norm<real_t>(rij);
        real_t dR = Rij - r0;

        if ( dR > 0.0 ) {
            energy_t energy{0.5 * fc * dR * dR};
            auto derHPdR = fc * dR;
            dist_vect_t unitVector = (Rij > eps ? (rij / Rij) : util::randomUnit());  // Should never occur.
            force_t f{};
            for (std::size_t k = 0; k != 3; ++k) {
                f[k] = -derHPdR * unitVector[k];
            }
            logger.trace("Leaving");
            return std::move(std::make_pair(energy, f));
        } else {
            logger.trace("Leaving.");
            return std::move(std::make_pair(energy_t{0.0}, force_t{0.0, 0.0, 0.0}));
        }
    }
}

