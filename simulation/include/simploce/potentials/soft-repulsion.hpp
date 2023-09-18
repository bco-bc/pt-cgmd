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
     * Standard conservative soft repulsion potential commonly used in dissipative particle dynamics. The force is
     * of the form A * (1.0 - r/r_c) * u where r <= r_c is the distance between two DPD particles and r_c is a cutoff
     * distance. The factor A represents the maximal strength of the interaction
     * and u is a unit vector.
     */
    class SoftRepulsion : public pair_potential {
    public:

        /**
         * Constructor.
         * @param forceField Force field.
         * @param bc Boundary conditions.
         * @param cutoffSR Cutoff distance for short range interactions.
         */
        SoftRepulsion(ff_ptr_t forceField, bc_ptr_t bc, dist_t cutoffSR);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

        /**
         * Returns potential energy and force on particle i.
         * @param rij Distance vector, r_ij = r_i - r_j, from j to i.
         * @param Rij Distance.
         * @param Rij2 Rij * Rij
         * @param Aij Maximum strength of interaction.
         * @param cutoffSR Cutoff distance for short ranged (SR) interactions.
         * @return Energy and force.
         */
        std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& rij,
                                                    real_t Rij,
                                                    real_t Rij2,
                                                    real_t Aij) const;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
        dist_t cutoffSR_;
    };
}

#endif //SIMULATION_SOFT_REPULSION_HPP
