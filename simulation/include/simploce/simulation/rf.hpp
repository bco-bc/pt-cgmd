/**
 * Author: André H. Juffer.
 * Created on 22/11/2021, 14:15
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_RF_HPP
#define SIMULATION_RF_HPP

#include "pair-potential.hpp"

namespace simploce {

    /**
     * Electrostatic interaction, reaction field field approach. This follows
     * Riniker et all, J. Chem. Phys., 134, 084110 (2011), specifically Eqs (2) - (5) without the
     * Lennard-Jones part. The last (constant) term of Eq (5) is excluded as well.
     * @see <a href="https://dx.doi.org/10.1063/1.3553378">
     * Riniker and van Gunsteren, J. Chem. Phys., 134, 084110.2011.
     * </a>
     */
    class RF : public pair_potential {
    public:

        /**
         * Constructor. All arguments are required.
         * @param kappa Inverse Debye screening length.
         * @param forceField Force field.
         * @param box Simulation box.
         * @param bc Boundary condition.
         */
        RF(real_t kappa, ff_ptr_t forceField, box_ptr_t box, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        friend class LJ_RF;

        std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& rij,
                                                    real_t Rij,
                                                    real_t Rij2,
                                                    charge_t q1,
                                                    charge_t q2,
                                                    real_t eps_inside_rc,
                                                    real_t eps_outside_rc);

        real_t compute_C_rf_(const distance_t& rc, real_t eps_inside_rc, real_t eps_outside_rc);

        real_t kappa_;
        ff_ptr_t forceField_;
        box_ptr_t box_;
        bc_ptr_t bc_;

    };
}

#endif //SIMULATION_RF_HPP
