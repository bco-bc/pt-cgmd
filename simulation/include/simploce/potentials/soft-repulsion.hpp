/*
 * Author: André H. Juffer.
 * Created on 23/05/2022, 17:23.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_SOFT_REPULSION_HPP
#define SIMULATION_SOFT_REPULSION_HPP

#include "simploce/potentials/pair-potential.hpp"

namespace simploce {

    /**
     * Standard conservative soft repulsion potential commonly used in dissipative particle dynamics.
     */
    class SoftRepulsion : public pair_potential {
    public:

        /**
         * Constructor.
         * @param forceField Force field.
         * @param bc Boundary conditions.
         * @param cutoff Cutoff distance.
         */
        SoftRepulsion(ff_ptr_t forceField, bc_ptr_t bc, dist_t cutoff);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

        /**
         * Returns potential energy and force on particle i.
         * @param r_ij Distance vector, r_ij = r_i - r_j, from j to i.
         * @param uv_ij Unit vector from j to i.
         * @param R_ij Distance.
         * @param a_ij Maximum strength of interaction.
         * @param cutoff Cutoff distance.
         * @return Energy and force.
         */
        static std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& r_ij,
                                                           const dist_vect_t& uv_ij,
                                                           real_t R_ij,
                                                           real_t a_ij,
                                                           const dist_t& cutoff);

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
        dist_t cutoff_;
    };
}

#endif //SIMULATION_SOFT_REPULSION_HPP
