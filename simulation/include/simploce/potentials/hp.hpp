/*
 * Author: André H. Juffer.
 * Created on 13/11/2021, 17:23.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_HP_HPP
#define SIMULATION_HP_HPP

#include "simploce/potentials/pair-potential.hpp"

namespace simploce {

    /**
     * Harmonic bonded potential, U(r) = 0.5 * k * (r - r0)^2, where r is a distance, k is a
     * force constant (fc) and r0 is the equilibrium distance.
     */
    class HP : public pair_potential {
    public:

        HP(ff_ptr_t forceField, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &pi, const p_ptr_t &pj) override;

        /**
         * Returns potential energy and force on particle i.
         * @param r_ij Distance vector, r_ij = r_i - r_j, from j to i.
         * @param uv_ij Unit vector from j to i.
         * @param R_ij Distance.
         * @param r0 Equilibrium distance.
         * @param fc Force constant.
         * @return Energy and force.
         */
        static std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t &r_ij,
                                                           const dist_vect_t& uv_ij,
                                                           real_t R_ij,
                                                           real_t r0,
                                                           real_t fc);

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
    };

}

#endif //SIMULATION_HP_HPP
