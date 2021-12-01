/*
 * Author: André H. Juffer.
 * Created on 23/11/2021, 12:31.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_SC_HPP
#define SIMULATION_SC_HPP

#include "simploce/potentials/pair-potential.hpp"

namespace simploce {

    /**
     * Screened Coulomb, i.e. 1.0/(4 * pi * eps_0 * eps_r * r.
     */
    class SC : public pair_potential {
    public:

        SC(ff_ptr_t forceField, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        friend class HS_SC;

        static std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& rij,
                                                           const real_t& Rij,
                                                           const real_t& Rij2,
                                                           const charge_t& q1,
                                                           const charge_t& q2,
                                                           real_t eps_inside_rc);

        ff_ptr_t forceField_;
        bc_ptr_t bc_;

    };
}

#endif //SIMULATION_SC_HPP
