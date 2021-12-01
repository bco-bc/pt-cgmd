/*
 * Author: André H. Juffer.
 * Created on 22/11/2021, 21:21.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_SF_HPP
#define SIMULATION_SF_HPP

#include "simploce/potentials/pair-potential.hpp"

namespace simploce {

    /**
     * Shifted force Coulomb interaction. The Coulomb interaction is calculated according to the
     * shifted force (SF) method of Levitt, M. et al, Comput. Phys. Commun. 1995, 91, 215−231.
     */
    class SF : public pair_potential {
    public:

        /**
         * Constructor. All argument are required.
         * @param forceField Force field.
         * @param box Simulation box.
         * @param bc Boundary condition.
         */
        SF(ff_ptr_t forceField, box_ptr_t box, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        friend class HS_SF;
        friend class LJ_SF;

        std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& rij,
                                                    real_t Rij,
                                                    real_t Rij2,
                                                    const charge_t& q1,
                                                    const charge_t& q2,
                                                    real_t eps_inside_rc);

        ff_ptr_t forceField_;
        box_ptr_t box_;
        bc_ptr_t bc_;

    };
}

#endif //SIMULATION_SF_HPP
