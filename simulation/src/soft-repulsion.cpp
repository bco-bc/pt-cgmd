/*
 * Author: André H. Juffer.
 * Created on 23/05/2022, 17:23.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/soft-repulsion.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <stdexcept>
#include <limits>

namespace simploce {

    SoftRepulsion::SoftRepulsion(ff_ptr_t forceField, bc_ptr_t bc, dist_t cutoffSR) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)}, cutoffSR_{cutoffSR} {
        util::Logger logger{"simploce::SoftRepulsion::SoftRepulsion(...)"};
        if (!forceField_) {
            logAndThrow(logger, "Missing force field.");
        }
        if (!bc_) {
            logAndThrow(logger, "Missing boundary conditions.");
        }
        if (cutoffSR_() <= 0.0) {
            logAndThrow(logger, "A cutoff distance must be > 0.0.");
        }
    }

    std::pair<energy_t, force_t>
    SoftRepulsion::operator()(const p_ptr_t &p1, const p_ptr_t &p2) {
        static util::Logger logger("simploce::SoftRepulsion::operator()");
        logger.trace("Entering");

        static auto eps = real_t(std::numeric_limits<float>::min());

        // Get interaction parameters.
        auto Aij = forceField_->softRepulsion(p1->spec(), p2->spec());

        // Current positions.
        const auto &ri = p1->position();
        const auto &rj = p2->position();

        // Apply boundary condition.
        auto rij = bc_->apply(ri, rj);
        auto Rij = norm<real_t>(rij);  // Should not be zero.
        auto Rij2 = Rij * Rij;

        // Forces and energy.
        auto pair = forceAndEnergy(rij, Rij, Rij2, Aij);

        // Done.
        logger.trace("Leaving");
        return std::move(pair);
    }

    std::pair<energy_t, force_t>
    SoftRepulsion::forceAndEnergy(const dist_vect_t &rij,
                                  real_t Rij,
                                  real_t Rij2,
                                  real_t Aij) const {
        static util::Logger logger("simploce::SoftRepulsion::forceAndEnergy()");
        logger.trace("Entering.");

        static auto eps = real_t(std::numeric_limits<float>::min());

        // Within cutoff distance?
        if (Rij >= cutoffSR_()) {
            return std::make_pair(energy_t{0.0}, force_t{0.0, 0.0, 0.0});
        } else if (Rij <= eps ) {
            logger.warn(std::to_string(Rij) + ": Zero distance encountered between two particles.");
        }

        // Interaction energy and force on particle i. Rij should never be zero.
        auto unitVector = (Rij > eps ? (rij / Rij) : util::randomUnit());
        energy_t energy{-Aij * (Rij - Rij2 / (2.0 * cutoffSR_())) + 0.5 * Aij * cutoffSR_()};
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = Aij * (1.0 - Rij / cutoffSR_ ()) * unitVector[k];
        }

        logger.trace("Leaving.");
        return std::move(std::make_pair(energy, f));
    }
}