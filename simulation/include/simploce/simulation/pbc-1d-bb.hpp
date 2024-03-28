/*
 * Author: André H. Juffer.
 * Created on 22/12/2021, 19:08.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_PBC_1D_BB_HPP
#define SIMULATION_PBC_1D_BB_HPP

#include "bc-impl.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {

    /**
     * Periodic boundary condition (PBC) in one dimension (1D) with the nearest image approximation,
     * as well as bounce-back (BB) reflection for velocities, that is a particle's velocity
     * components are reversed (same magnitude but opposite direction) if it crossed the
     * box's boundary.
     * @see Revenga et al, International Journal of Modern Physics C, v. 9, p. 1319, 1998.,
     * Revenga et al, Computer Physics Communications, v. 121-122, p. 309 - 311, 1999.
     */
    class PBC_1D_BB : public boundary_condition_impl {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param direction Direction for periodicity is to be applied.
         */
        explicit PBC_1D_BB(box_ptr_t box,
                           Direction direction = Direction::Z);

        dist_vect_t
        apply(const position_t& ri,
              const position_t& rj) const override;

        position_t
        placeInside(const position_t& r_out) const override;

        /**
         * No effect.
         * @return r
         */
        position_t
        apply(const position_t& r) const override;

        /**
         * Applies bounce-back.
         * @param v Current velocity.
         * @param r Current position.
         * @return Bounced-back velocity.
         */
        velocity_t
        apply(const velocity_t &v,
              const position_t& r) const override;

        /**
         * Applies bounce-back to all velocities of all particles in the group, if the group
         * has crossed the box boundaries.
         * @param particleGroup Particle group.
         */
        void
        applyToVelocities(const pg_ptr_t& particleGroup) const override;

    private:

        box_ptr_t box_;
        Direction direction_;
    };

}

#endif //SIMULATION_PBC_1D_BB_HPP
