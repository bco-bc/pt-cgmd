/*
 * Author: André H. Juffer.
 * Created on 08/08/2023
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */


#ifndef SIMULATION_GAUSS_SF_SR_HPP
#define SIMULATION_GAUSS_SF_SR_HPP

#include "pair-potential.hpp"
#include "gauss-sf.hpp"
#include "soft-repulsion.hpp"

namespace simploce {

    /**
     * Interaction between overlapping Gaussian charge densities plus soft-repulsion.
     * This potential is for mesoscopic simulations.
     */
    class GaussianSF_SoftRepulsion : public pair_potential {
    public:

        /**
         * Constructor.
         * @param forceField Force field.
         * @param bc Boundary conditions.
         * @param cutoffs Cutoff distances.
         * @param gaussianSF Pair potential for overlapping Gaussian charge densities.
         * @param softRepulsion Pair potential for the soft repulsive interaction.
         *
         */
        GaussianSF_SoftRepulsion(ff_ptr_t forceField,
                                 bc_ptr_t bc,
                                 std::shared_ptr<GaussianSF> gaussianSF,
                                 std::shared_ptr<SoftRepulsion> softRepulsion);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &pi, const p_ptr_t &pj) override;

        std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& rij,
                                                    real_t Rij,
                                                    real_t Rij2,
                                                    real_t Aij,
                                                    const charge_t& q_i,
                                                    const charge_t& q_j,
                                                    real_t sigma_i,
                                                    real_t sigma_j);

    private:
        ff_ptr_t forceField_;
        bc_ptr_t bc_;

        std::shared_ptr<GaussianSF> gaussianSF_;
        std::shared_ptr<SoftRepulsion> softRepulsion_;
    };
}

#endif //SIMULATION_GAUSS_SF_SR_HPP
