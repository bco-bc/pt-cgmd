/*
 * Author: André H. Juffer.
 * Created on 07/08/2023
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/solid-sphere-dsf.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"

namespace simploce {

    SolidSphere_DSF::SolidSphere_DSF(simploce::dist_t cutoff,
                                     simploce::ff_ptr_t forceField,
                                     simploce::box_ptr_t box,
                                     simploce::bc_ptr_t bc,
                                     simploce::dist_t radius,
                                     bool mesoscopic) :
        cutoff_{cutoff}, forceField_{std::move(forceField)}, box_{std::move(box)}, bc_{std::move(bc)},
        radius_{radius}, mesoscopic_{mesoscopic} {
    }

    std::pair<energy_t, force_t>
    SolidSphere_DSF::operator()(const simploce::p_ptr_t &p1, const simploce::p_ptr_t &p2) {
         // Current positions.
        const auto &r1 = p1->position();
        const auto &r2 = p2->position();

        // Charges
        auto q1 = p1->charge();
        auto q2 = p2->charge();

        // Apply boundary conditions.
        dist_vect_t rij = bc_->apply(r1, r2);
        auto Rij = norm<real_t>(rij);
        real_t Rij2 = Rij * Rij;

        // Potential/interaction energy.
        return std::move(this->forceAndEnergy(rij, Rij, Rij2, q1, q2));
    }

    std::pair<energy_t, force_t>
    SolidSphere_DSF::forceAndEnergy(const simploce::dist_vect_t &rij,
                                    simploce::real_t Rij,
                                    simploce::real_t Rij2,
                                    const simploce::charge_t &q1,
                                    const simploce::charge_t &q2) {
        static util::Logger logger("simploce::SoftRepulsion::forceAndEnergy()");
        logger.trace("Entering.");

        static auto rc2 = cutoff_() * cutoff_();
        static auto radius2 = radius_() * radius_();
        static auto radius3 = radius2 * radius_();
        static auto eps = real_t(std::numeric_limits<float>::min());

        if (Rij2 > rc2) {
            return std::move(std::make_pair(energy_t{}, force_t{}));
        }

        real_t t1 = q1() * q2() / (units::mu<real_t>::FOUR_PI_E0);
        dist_vect_t uv = Rij > eps ? rij/Rij : rij;
        energy_t energy{};
        force_t f{};
        if (Rij2 > radius2) {
            // Outside solid sphere.
            real_t dElecdR = t1 * (-1.0 / Rij2 + 1.0 / rc2);
            energy = t1 * (1.0 / Rij - 1.0 / cutoff_() + (Rij - cutoff_()) / rc2);
            for (std::size_t k = 0; k != 3; ++k) {
                f[k] = -dElecdR * uv[k];
            }
        } else {
            // Inside solid sphere.
            auto pot_r = 0.5 * (3.0 - Rij2 / radius2) / radius_();
            real_t dElecdR = t1 * (-Rij / radius3 + 1.0 / rc2);
            energy = t1 * (pot_r - 1.0 / cutoff_() + (Rij - cutoff_()) / rc2);
            for (std::size_t k = 0; k != 3; ++k) {
                f[k] = -dElecdR * uv[k];
            }
        }

        logger.trace("Leaving.");
        return std::move(std::make_pair(energy, f));
    }
}
