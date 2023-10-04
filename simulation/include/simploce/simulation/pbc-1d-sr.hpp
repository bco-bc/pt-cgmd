/*
 * Author: André H. Juffer.
 * Created on 26/06/2023.
 *
 * Copyright (c) Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_PBC_1D_SR_HPP
#define SIMULATION_PBC_1D_SR_HPP

#include "bc.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {

    /**
     * Periodic boundary conditions (PBC) in one dimension (1D) with the nearest image approximation, as well as
     * specular reflection (SR) for velocities. That is, the parallel component of the particle's velocity is retained,
     * but the normal component (perpendicular to the direction) is reversed (same magnitude but opposite value),
     * if the particle crossed the box's boundary.
     * @see Revenga et al, International Journal of Modern Physics C, v. 9, p. 1319, 1998.,
     * Revenga et al, Computer Physics Communications, v. 121-122, p. 309 - 311, 1999.
     */
    class PBC_1D_SR: public boundary_condition {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param direction Direction where periodicity is to be applied.
         */
        explicit PBC_1D_SR(box_ptr_t box,
                           Direction direction = Direction::Z);

        dist_vect_t apply(const position_t& ri, const position_t& rj) const override;

        position_t placeInside(const position_t& r_out) const override;

        /**
         * Applies specular reflection.
         * @param v Current velocity.
         * @param r Current position.
         * @return Reflected velocity.
         */
        velocity_t apply(const velocity_t &v,
                         const position_t& r) const override;

        /**
         * Applies specular reflection to all velocities of all particles in the group, if the group
         * has crossed the box boundaries.
         * @param particleGroup Particle group.
         */
        void applyToVelocities(const pg_ptr_t& particleGroup) const override;

    private:

        box_ptr_t box_;
        Direction direction_;
    };
}

#endif //SIMULATION_PBC_1D_SR_HPP
