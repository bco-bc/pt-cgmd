/*
 * Author: André H. Juffer.
 * Created on 08/08/2023
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/gauss-sf.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/math-constants.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <utility>

namespace simploce {

    /**
     * Returns factor S for erf(x) where x = S*r
     * @param sigma_i Width density.
     * @param sigma_j Width density.
     * @return Value of S.
     */
    static real_t
    S_(real_t sigma_i,
       real_t sigma_j) {
        return std::sqrt(1.0 / (sigma_i * sigma_i + sigma_j * sigma_j));
    }

    /**
     * Returns U(r) and dU(r)/dr at cutoff distance of the Gaussian interaction potential U(r).
     * @param q_i Charge value
     * @param q_j Charge value
     * @param S Factor for error function.
     * @param cutoffLR Cutoff distance.
     * @param mesoscopic For a mesoscopic interaction?
     * @return U(r), dU(r)/du at cutoff distance.
     */
    static std::pair<real_t, real_t>
    atCutoff(const charge_t& q_i,
             const charge_t& q_j,
             real_t S,
             dist_t cutoffLR,
             bool mesoscopic) {
        static auto four_pi_e0 = units::mu<real_t>::FOUR_PI_E0;
        static auto four_pi = 4.0 * math::constants<real_t>::PI;
        static auto sqrt_pi = std::sqrt(math::constants<real_t>::PI);

        auto r = cutoffLR();
        auto x = S * r;
        auto erf_x = std::erf(x);
        auto t = erf_x * q_i() * q_j() / cutoffLR();
        auto u = mesoscopic ? t / four_pi : t / four_pi_e0;

        auto e = std::exp(-x * x);
        t = (2.0/sqrt_pi * e * x - erf_x) * q_i() * q_j() / (r * r);
        auto dUdr = mesoscopic ? t / four_pi : t / four_pi_e0;

        return std::move(std::make_pair(u, dUdr));
    }

    GaussianSF::GaussianSF(simploce::ff_ptr_t forceField, simploce::bc_ptr_t bc, dist_t cutoffLR, bool mesoscopic) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)}, cutoffLR_{cutoffLR}, mesoscopic_{mesoscopic} {
        util::Logger logger("simploce::GaussianSF::GaussianSF(...)");
        if (!forceField_) {
            logAndThrow(logger, "Missing force field.");
        }
        if (!bc_) {
            logAndThrow(logger, "Missing boundary conditions.");
        }
        if (cutoffLR_() <= 0) {
            logAndThrow(logger, "A cutoff distance must be > 0.0.");
        }
    }

    std::pair<energy_t, force_t>
    GaussianSF::operator()(const simploce::p_ptr_t &pi, const simploce::p_ptr_t &pj) {
        static util::Logger logger("simploce::GaussianSF::operator()(...)");
        logger.trace("Entering");

        // Get interaction parameters.
        auto params = forceField_->gaussianChargeDensity(pi->spec(), pj->spec());
        auto sigma_i = params.first;
        auto sigma_j = params.second;
        auto q_i = pi->charge();
        auto q_j = pj->charge();

        // Current positions.
        const auto &ri = pi->position();
        const auto &rj = pj->position();

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(ri, rj);
        auto Rij = norm<real_t>(rij);
        real_t Rij2 = Rij * Rij;

        // Forces and energy.
        auto pair = std::move(this->forceAndEnergy(rij,
                                                  Rij,
                                                  Rij2,
                                                  q_i,
                                                  q_j,
                                                  sigma_i,
                                                  sigma_j));

        // Done.
        logger.trace("Leaving");
        return std::move(pair);
    }

    std::pair<energy_t, force_t>
    GaussianSF::forceAndEnergy(const dist_vect_t &rij,
                               real_t Rij,
                               real_t Rij2,
                               const charge_t& q_i,
                               const charge_t& q_j,
                               real_t sigma_i,
                               real_t sigma_j) {
        static util::Logger logger("simploce::GaussianSF::forceAndEnergy(...)");
        logger.trace("Entering");

        static auto four_pi_e0 = units::mu<real_t>::FOUR_PI_E0;
        static auto four_pi = 4.0 * math::constants<real_t>::PI;
        static auto sqrt_pi = std::sqrt(math::constants<real_t>::PI);
        static auto eps = real_t(std::numeric_limits<float>::min());

        // Within cutoff distance?
        if (Rij >= cutoffLR_()) {
            return std::make_pair(energy_t{0.0}, force_t{0.0, 0.0, 0.0});
        } else if (Rij <= eps ) {
            logger.warn(std::to_string(Rij) + ": Zero distance encountered between two particles.");
        }

        // Interaction energy and force on particle i. Rij should never be zero.
        auto S = S_(sigma_i, sigma_j);
        auto pair = atCutoff(q_i, q_j, S, cutoffLR_, mesoscopic_);
        auto UAtCutoff = pair.first;
        auto dUdrAtCutoff = pair.second;

        auto x = S * Rij;
        auto erf_x = std::erf(x) ;
        auto t = erf_x * q_i() * q_j() / Rij;
        t  = mesoscopic_ ? t / four_pi : t / four_pi_e0;
        auto energy = t - UAtCutoff - dUdrAtCutoff * (Rij - cutoffLR_());  // Shifted force

        auto e = std::exp(-x * x);
        t = (2.0 / sqrt_pi * e * x - erf_x) * q_i() * q_j() / Rij2;
        auto dUdr = mesoscopic_ ? t / four_pi : t / four_pi_e0;
        dist_vect_t unitVector = (Rij > eps ? (rij / Rij) : util::randomUnit());
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = (-dUdr + dUdrAtCutoff) * unitVector[k];                 // Shifted force
        }

        // Done
        logger.trace("Leaving");
        return std::move(std::make_pair(energy, f));
    }
}
