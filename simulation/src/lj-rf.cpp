/**
 * Author: André H. Juffer.
 * Created on 15/11/2021, 14:20
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/lj-rf.hpp"
#include "simploce/simulation/force-field.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"

namespace simploce {

    LJ_RF::LJ_RF(real_t kappa, ff_ptr_t forceField, box_ptr_t box, bc_ptr_t bc) :
        kappa_(kappa), forceField_(std::move(forceField)), box_{std::move(box)},
        bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    LJ_RF::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        // Initialize.
        static const distance_t rc = properties::cutoffDistance(box_);
        static const real_t rc2 = rc() * rc();
        static const distance_t rf = rc;
        static const real_t rf3 = rf() * rf() * rf();
        static const real_t C_rf = compute_C_rf_(rc);
        static const real_t eps_cs = forceField_->relativePermittivityInsideCutoff();

        // Get interaction parameters.
        auto params = forceField_->lennardJonesReactionFieldParameters(p1->spec(), p2->spec());
        auto C12 = std::get<0>(params);
        auto C6 = std::get<1>(params);

        // Current positions.
        const auto &r1 = p1->position();
        const auto &r2 = p2->position();

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(r1, r2);
        auto Rij = norm<real_t>(rij);

        // Charges
        auto q1 = p1->charge();
        auto q2 = p2->charge();

        // Potential/interaction energy, kJ/mol.
        real_t Rij2 = Rij * Rij;
        real_t Rij6 = Rij2 * Rij2 * Rij2;
        real_t Rij12 = Rij6 * Rij6;
        real_t t1 = C12 / Rij12;
        real_t t2 = C6 / Rij6;
        real_t c1 = q1() * q2() / (units::mu<real_t>::FOUR_PI_E0 * eps_cs);
        real_t c2 = -c1 * (0.5 * C_rf * Rij2 / rf3 + (1.0 - 0.5 * C_rf) / rf());
        energy_t energy{t1 - t2 + c1 / Rij + c2};  // LJ

         // Force on particle #1, kJ/(mol nm).
        dist_vect_t unitVector = rij / Rij;
        real_t dLJdR = -6.0 * ( 2.0 * t1 - t2 ) / Rij;
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = -dLJdR * unitVector[k];
        }

        // Done.
        return std::move(std::make_pair(energy, f));
    }

    real_t LJ_RF::compute_C_rf_(const distance_t& rc) {
        real_t kappa_rc = kappa_ * rc();
        real_t kappa_rc_2 = kappa_rc * kappa_rc;
        real_t eps_rf = (1.0 + kappa_rc_2 / (2.0 * (kappa_rc + 1.0)) * forceField_->relativePermittivityBeyondCutoff());
        real_t c1 = 2.0 * forceField_->relativePermittivityInsideCutoff() - 2.0 * eps_rf;
        real_t c2 = forceField_->relativePermittivityInsideCutoff() + 2.0 * eps_rf;
        real_t c3 = 1.0 + kappa_rc;
        return (c1 * c3 - eps_rf * kappa_rc_2) / (c2 * c3 + eps_rf * kappa_rc_2);
    }

}