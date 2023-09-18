/*
 * Author: André H. Juffer.
 * Created on 08/08/2023
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/gauss-sf-sr.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/util/logger.hpp"

namespace simploce {

    GaussianSF_SoftRepulsion::GaussianSF_SoftRepulsion(ff_ptr_t forceField,
                                                       bc_ptr_t bc,
                                                       std::shared_ptr<GaussianSF> gaussianSF,
                                                       std::shared_ptr<SoftRepulsion> softRepulsion) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)},
        gaussianSF_{std::move(gaussianSF)}, softRepulsion_{std::move(softRepulsion)} {
        util::Logger logger{"simploce::GaussianSF_SoftRepulsion::GaussianSF_SoftRepulsion()"};
        if (!forceField_) {
            logAndThrow(logger, "Missing force field.");
        }
        if (!bc_) {
            logAndThrow(logger, "Missing boundary conditions.");
        }
        if (!gaussianSF_) {
            logAndThrow(logger, "Potential for interacting Gaussian densities is missing.");
        }
        if (!softRepulsion_) {
            logAndThrow(logger, "Potential for soft repulsion is missing.");
        }
    }

    std::pair<energy_t, force_t>
    GaussianSF_SoftRepulsion::operator()(const p_ptr_t &pi,
                                         const p_ptr_t &pj) {
        static util::Logger logger("simploce::GaussianSF_SoftRepulsion::operator()()");
        logger.trace("Entering.");

        auto params =
                forceField_->gaussianChargeDensitySoftRepulsion(pi->spec(), pj->spec());
        real_t sigma_i = std::get<0>(params);
        real_t sigma_j = std::get<1>(params);
        real_t Aij = std::get<2>(params);

        // Current positions.
        const auto &ri = pi->position();
        const auto &rj = pj->position();

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(ri, rj);
        auto Rij = norm<real_t>(rij);
        real_t Rij2 = Rij * Rij;
        auto pair = forceAndEnergy(rij, Rij, Rij2, Aij, pi->charge(), pj->charge(), sigma_i, sigma_j);

        // Done
        logger.trace("Leaving.");
        return std::move(pair);
    }

    std::pair<energy_t, force_t>
    GaussianSF_SoftRepulsion::forceAndEnergy(const dist_vect_t &rij,
                                             real_t Rij,
                                             real_t Rij2,
                                             real_t Aij,
                                             const charge_t &q_i,
                                             const charge_t &q_j,
                                             real_t sigma_i,
                                             real_t sigma_j) {
        static util::Logger logger{"simploce::GaussianSF_SoftRepulsion::forceAndEnergy()"};
        logger.trace("Entering.");

        // Combine energies and forces from the other potentials.
        auto pairGaussianSF =
                gaussianSF_->forceAndEnergy(rij, Rij, Rij2, q_i, q_j, sigma_i, sigma_j);
        auto pairSR = softRepulsion_->forceAndEnergy(rij, Rij, Rij2, Aij);
        energy_t energy{pairGaussianSF.first + pairSR.first};
        force_t f{pairGaussianSF.second + pairSR.second};

        // Done.
        logger.trace("Leaving.");
        return std::move(std::make_pair(energy, f));
    }
}
