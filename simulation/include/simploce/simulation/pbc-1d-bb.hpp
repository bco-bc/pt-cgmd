/*
 * Author: André H. Juffer.
 * Created on 22/12/2021, 19:08.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_PBC_1D_BB_HPP
#define SIMULATION_PBC_1D_BB_HPP

#include "bc.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {

    /**
     * Periodic boundary condition (PBC) in one dimension (1D) with the nearest image approximation, as well as
     * bounce-back (BB) reflection for velocities, that is a particle's velocity components are reversed
     * (same magnitude but opposite direction) if it crossed the boundary.
     */
    class PBC_1D_BB : public boundary_condition {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param direction Direction where periodicity is to be applied.
         */
        explicit PBC_1D_BB(box_ptr_t box,
                           Direction direction = Direction::Z);

        dist_vect_t apply(const position_t& ri, const position_t& rj) const override;

        position_t placeInside(const position_t& r_out) const override;

        /**
         * Applies bounce-back.
         * @param v Current velocity.
         * @param r Current position.
         * @return Bounced-back velocity.
         */
        velocity_t apply(const velocity_t &v,
                         const position_t& r) const override;

    private:

        box_ptr_t box_;
        Direction direction_;
    };

}

#endif //SIMULATION_PBC_1D_BB_HPP
