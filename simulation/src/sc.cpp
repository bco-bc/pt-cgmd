/*
 * Author: André H. Juffer.
 * Created on 23/11/2021, 12:39.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/sc.hpp"
#include "simploce/simulation/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"

#include <utility>

namespace simploce {

    SC::SC(ff_ptr_t forceField, bc_ptr_t bc) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    SC::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
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

        // Interaction parameters.
        auto eps_inside_rc = forceField_->screenedCoulomb(p1->spec(), p2->spec());

        // Potential/interaction energy, kJ/mol.
        return std::move(SC::forceAndEnergy(rij, Rij, Rij2, q1, q2, eps_inside_rc));
    }

    std::pair<energy_t, force_t>
    SC::forceAndEnergy(const dist_vect_t& rij,
                       const real_t& Rij,
                       const real_t& Rij2,
                       const charge_t& q1,
                       const charge_t& q2,
                       real_t eps_inside_rc) {
        auto c1 = 1.0 / (units::mu<real_t>::FOUR_PI_E0 * eps_inside_rc);
        energy_t energy{c1 * q1() * q2() / Rij};
        auto dSFdr = -c1 / Rij;
        force_t f{};
        dist_vect_t unitVector = rij / Rij;
        for (int k = 0; k !=3; ++k) {
            f[k] = -dSFdr * unitVector[k];
        }

        return std::move(std::make_pair(energy, f));
    }
}

