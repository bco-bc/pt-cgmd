/**
 * Author: André H. Juffer.
 * Created on 15/11/2021, 14:20
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_HALVE_ATTRACTIVE_QP_HPP
#define SIMULATION_HALVE_ATTRACTIVE_QP_HPP

#include "simploce/potentials/pair-potential.hpp"
#include "simploce/types/s-types.hpp"

namespace simploce {

    /**
     * Halve attractive quartic potential, U(r) = 0.5 * k * (r - r0)^4, where r is a distance,
     * k is a force constant (fc) and r0 is the equilibrium distance. Halve attractive implies
     * U(r) = 0 for r <= r0.
     * @see <a href="https://dx.doi.org/10.1063/1.3553378">
     * Riniker and van Gunsteren, J. Chem. Phys., 134, 084110.2011.
     * </a>
     */
    class HalveAttractiveQP : public pair_potential {
    public:

        HalveAttractiveQP(ff_ptr_t forceField, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
    };


}

#endif //SIMULATION_HALVE_ATTRACTIVE_QP_HPP
