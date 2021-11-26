/**
 * Author: André H. Juffer.
 * Created on 15/11/2021, 14:20
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_LJ_RF_HPP
#define SIMULATION_LJ_RF_HPP

#include "pair-potential.hpp"
#include "rf.hpp"
#include "s-types.hpp"
#include <utility>

namespace simploce {

    /**
     * Lennard-Jones interaction plus reaction field (RF) electrostatics. This follows
     * Riniker et all, J. Chem. Phys., 134, 084110 (2011), specifically Eqs (2) - (5).
     * The last constant term of Eq (5) is excluded.
     * @see <a href="https://dx.doi.org/10.1063/1.3553378">
     * Riniker and van Gunsteren, J. Chem. Phys., 134, 084110.2011.
     * </a>
     */
    class LJ_RF : public pair_potential {
    public:

        /**
         * Constructor. All arguments are required.
         * @param kappa Inverse Debye screening length.
         * @param forceField Force field.
         * @param box Simulation box.
         * @param bc Boundary condition.
         */
        LJ_RF(real_t kappa, ff_ptr_t forceField, box_ptr_t box, bc_ptr_t bc, rf_ptr_t  rf);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
        rf_ptr_t rf_;

    };

}

#endif //SIMULATION_LJ_RF_HPP
