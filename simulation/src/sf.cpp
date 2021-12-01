/*
 * Author: André H. Juffer.
 * Created on 22/11/2021, 21:25.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/sf.hpp"

#include <utility>
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/units/units-mu.hpp"

namespace simploce {

    SF::SF(ff_ptr_t forceField, box_ptr_t box, bc_ptr_t bc) :
        forceField_{std::move(forceField)}, box_{std::move(box)}, bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    SF::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
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
        real_t eps_inside_rc = forceField_->shiftedForce(p1->spec(), p2->spec());

        // Potential/interaction energy, kJ/mol.
        return std::move(this->forceAndEnergy(rij, Rij, Rij2, q1, q2, eps_inside_rc));
    }

    std::pair<energy_t, force_t>
    SF::forceAndEnergy(const dist_vect_t& rij,
                       real_t Rij,
                       real_t Rij2,
                       const charge_t& q1,
                       const charge_t& q2,
                       real_t eps_inside_rc) {
        static auto rc = properties::cutoffDistance(box_);
        static auto rc2 = rc() * rc();

        real_t t1 = q1() * q2() / (units::mu<real_t>::FOUR_PI_E0 * eps_inside_rc);
        real_t elec = t1 * (1.0 / Rij - 1.0 / rc() + (Rij - rc()) / rc2);    // kJ/mol

        dist_vect_t uv = rij/Rij;
        real_t dElecdR = t1 * (-1.0 / Rij2 + 1.0 / rc2);
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = -dElecdR * uv[k];   // kJ/(mol nm)
        }

        return std::move(std::make_pair(elec, f));
    }
}

