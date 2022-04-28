/*
 * Author: André H. Juffer.
 * Created on 22/12/2021, 19:08.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_1D_PBC_HPP
#define SIMULATION_1D_PBC_HPP

#include "bc.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {

    /**
     * Periodic boundary condition in one dimension with the nearest image approximation.
     */
    class OneD_PeriodicBoundaryCondition : public boundary_condition {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param direction Direction where periodicity is to be applied.
         */
        explicit OneD_PeriodicBoundaryCondition(box_ptr_t box, Direction direction = Direction::Z);

        dist_vect_t apply(const position_t& ri, const position_t& rj) const override;

        position_t placeInside(const position_t& r_out) const override;

    private:

        box_ptr_t box_;
        Direction direction_;
    };

}

#endif //SIMULATION_1D_PBC_HPP
