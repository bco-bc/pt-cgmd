/**
 * Author: André H. Juffer.
 * Created on 22/11/2021, 14:21
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/rf.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <utility>

namespace simploce {

    RF::RF(real_t kappa, ff_ptr_t forceField, box_ptr_t box, bc_ptr_t bc) :
        kappa_(kappa), forceField_(std::move(forceField)), box_{std::move(box)},
        bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    RF::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
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
        auto params = forceField_->reactionField(p1->spec(), p2->spec());

        // Potential/interaction energy, kJ/mol.
        return std::move(this->forceAndEnergy(rij, Rij, Rij2, q1, q2, params.first, params.second));
    }

    std::pair<energy_t, force_t>
    RF::forceAndEnergy(const dist_vect_t& rij,
                       real_t Rij,
                       real_t Rij2,
                       charge_t q1,
                       charge_t q2,
                       real_t eps_inside_rc,
                       real_t eps_outside_rc) {
        static const dist_t rc = properties::cutoffDistance(box_);
        static const real_t rf = rc();
        static const real_t rf3 = rf * rf * rf;

        auto C_rf = compute_C_rf_(rc, eps_inside_rc, eps_outside_rc);
        auto factor = units::mu<real_t>::FOUR_PI_E0 * eps_inside_rc;

        real_t c1 = q1() * q2() / factor;
        real_t C = c1 / Rij;
        real_t RF = -c1 * (0.5 * C_rf * Rij2 / rf3 + (1.0 - 0.5 * C_rf) / rf);
        energy_t energy{C + RF};

         // Force on particle #1, kJ/(mol nm).
        dist_vect_t unitVector = rij / Rij;
        real_t dCRFdR = 0.0;
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = -dCRFdR * unitVector[k];
        }
        return std::move(std::pair<energy_t, force_t>{energy, f});
    }

    real_t
    RF::compute_C_rf_(const dist_t& rc, real_t eps_inside_rc, real_t eps_outside_rc) {
        util::Logger logger("simploce::RF::operator ()");
        real_t kappa_rc = kappa_ * rc();
        real_t kappa_rc_2 = kappa_rc * kappa_rc;

        // See Eq (33) of Christen et al., J. Comput. Chem. 26: 1719 –1751, 2005
        real_t eps_rf =
                (1.0 + kappa_rc_2 / (2.0 * (kappa_rc + 1.0))) * eps_outside_rc;

        real_t c1 = 2.0 * eps_inside_rc - 2.0 * eps_rf;
        real_t c2 = eps_inside_rc + 2.0 * eps_rf;
        real_t c3 = 1.0 + kappa_rc;
        real_t C_rf = (c1 * c3 - eps_rf * kappa_rc_2) / (c2 * c3 + eps_rf * kappa_rc_2);
        logger.debug("eps_cs: " + util::toString(eps_inside_rc));
        logger.debug("kappa: " + util::toString(kappa_));
        logger.debug("eps_rf: " + util::toString(eps_rf));
        logger.debug("C_rf: " + util::toString(C_rf));
        return C_rf;
    }

}

