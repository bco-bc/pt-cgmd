/*
 * Author: André H. Juffer.
 * Created on 08/08/2023
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */


#ifndef SIMULATION_GAUSS_SF_HPP
#define SIMULATION_GAUSS_SF_HPP

#include "simploce/potentials/pair-potential.hpp"

namespace simploce {

    /**
     * Interaction between Gaussian charge densities, where the "shifted force" method
     * is applied such that the interaction energy and associated forces are
     * exactly zero at distances equal to or larger than a given cutoff distance.
     * The interaction is "soft" at very short distances. The Gaussian interaction
     * energy U(r) is of the for U(r)=erf(Sr) * Q1 * Q2 / (4 * pi * eps0 * r),
     * where Q1 and Q2 are charge values of the charge densities, erf(x) is the
     * error function, pi is 3.14, eps0 is the vacuum permittivity, and r is distance.
     * The factor S is given from (1/(sigma1^2 + sigma2^2)^1/2, where sigma1(2) is
     * the width of the charge density of Q1(2).
     * @see Fennel and Gezelter, J. Chem. Phys. 124, 234104, 2006.
     */
    class GaussianSF : public pair_potential {
    public:

        /**
         * Constructor
         * @param forceField Force field.
         * @param bc Boundary condition.
         * @param cutoffLR Cutoff distance for long-ranged (LR) interaction.
         * @param mesoscopic If true, this potential is for mesoscopic simulations (e.g., DPD).
         */
        GaussianSF(ff_ptr_t forceField, bc_ptr_t bc, dist_t cutoffLR, bool mesoscopic=false);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &pi, const p_ptr_t &pj) override;

        std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& rij,
                                                    real_t Rij,
                                                    real_t Rij2,
                                                    const charge_t& q_i,
                                                    const charge_t& q_j,
                                                    real_t sigma_i,
                                                    real_t sigma_j);

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
        dist_t cutoffLR_;
        bool mesoscopic_;
    };
}

#endif //SIMULATION_GAUSS_SF_HPP
