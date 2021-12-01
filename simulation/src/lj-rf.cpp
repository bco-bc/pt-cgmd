/**
 * Author: André H. Juffer.
 * Created on 15/11/2021, 14:20
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/lj-rf.hpp"
#include "simploce/potentials/lj.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"
#include <utility>

namespace simploce {

    LJ_RF::LJ_RF(real_t kappa, ff_ptr_t forceField, box_ptr_t box, bc_ptr_t bc, rf_ptr_t  rf) :
        forceField_(std::move(forceField)), bc_{std::move(bc)}, rf_{std::move(rf)} {
    }

    // Self term in Eq (5) of Riniker et al. is excluded.
    std::pair<energy_t, force_t>
    LJ_RF::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        // Get interaction parameters.
        auto params = forceField_->lennardJonesReactionField(p1->spec(), p2->spec());
        auto C12 = std::get<0>(params);
        auto C6 = std::get<1>(params);
        auto eps_inside_rc = std::get<2>(params);
        auto eps_outside_rc = std::get<3>(params);

        // Current positions.
        const auto &r1 = p1->position();
        const auto &r2 = p2->position();

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(r1, r2);
        auto Rij = norm<real_t>(rij);
        real_t Rij2 = Rij * Rij;

        // Charges
        auto q1 = p1->charge();
        auto q2 = p2->charge();

        // Potential/interaction energy, kJ/mol.
        // Lennard-Jones
        auto lj = LJ::forceAndEnergy(rij, Rij, Rij2, C12, C6);
        // Coulomb + Reaction field.
        auto el = rf_->forceAndEnergy(rij, Rij, Rij2, q1, q2, eps_inside_rc, eps_outside_rc);
        energy_t energy{lj.first() + el.first()};

         // Force on particle #1, kJ/(mol nm).
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = lj.second[k] + el.second[k];
        }

        // Done.
        return std::move(std::make_pair(energy, f));
    }

}