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

    SoftRepulsion::SoftRepulsion(ff_ptr_t forceField, bc_ptr_t bc, dist_t cutoff) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)}, cutoff_{cutoff} {
        if (!forceField_) {
            throw std::domain_error("SoftRepulsion: Missing force field.");
        }
        if (!bc_) {
            throw std::domain_error("SoftRepulsion: Missing boundary conditions.");
        }
        if (cutoff_ <= 0.0) {
            throw std::domain_error("SoftRepulsion: Cutoff distance must be > 0.0.");
        }
    }

    std::pair<energy_t, force_t>
    SoftRepulsion::operator()(const p_ptr_t &p1, const p_ptr_t &p2) {
        static util::Logger logger("simploce::SoftRepulsion::operator()");
        logger.trace("Entering");

        static auto eps = real_t(std::numeric_limits<float>::min());
        logger.debug(p1->name() + ", " + p2->name() + ": Particle pair.");
        auto a_ij = forceField_->softRepulsion(p1->spec(), p2->spec());
        logger.debug(std::to_string(a_ij) + ": Maximum repulsion.");

        // Current positions.
        const auto &ri = p1->position();
        const auto &rj = p2->position();

        // Apply boundary condition.
        auto rij = bc_->apply(ri, rj);
        auto Rij = norm<real_t>(rij);  // Should not be zero.
        auto unitVector = (Rij > eps ? (rij / Rij) : util::randomUnit());

        logger.trace("Leaving.");
        return std::move(SoftRepulsion::forceAndEnergy(rij, unitVector, Rij, a_ij, cutoff_));
    }

    std::pair<energy_t, force_t>
    SoftRepulsion::forceAndEnergy(const dist_vect_t &r_ij,
                                  const dist_vect_t& uv_ij,
                                  real_t R_ij,
                                  real_t a_ij,
                                  const dist_t &cutoff) {
        static util::Logger logger("simploce::SoftRepulsion::forceAndEnergy()");
        logger.trace("Entering.");

        static auto eps = real_t(std::numeric_limits<float>::min());

        // Within cutoff distance?
        if (R_ij >= cutoff()) {
            return std::make_pair(energy_t{0.0}, force_t{0.0, 0.0, 0.0});
        } else if (R_ij <= eps ) {
            logger.warn(std::to_string(R_ij) + ": Zero distance encountered between two particles.");
        }
        // Interaction energy and force on particle i. R_ij should never be zero.
        energy_t energy{-a_ij * (R_ij - R_ij * R_ij / (2.0 * cutoff())) + 0.5 * a_ij * cutoff()};
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = a_ij * (1.0 - R_ij / cutoff()) * uv_ij[k];
        }

        logger.trace("Leaving.");
        return std::move(std::make_pair(energy, f));
    }
}