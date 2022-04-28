/*
 * Author: André H. Juffer.
 * Created on 22/12/2021, 19:14.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/1d-pbc.hpp"
#include "simploce/util/util.hpp"

namespace simploce {

    static int index_(const Direction& direction) {
        if ( direction == Direction::X) {
            return 0;
        } else if (direction == Direction::Y) {
            return 1;
        } else {
            return 2;
        }
    }

    OneD_PeriodicBoundaryCondition::OneD_PeriodicBoundaryCondition(box_ptr_t box, Direction direction) :
        box_{std::move(box)}, direction_{direction} {
    }

    dist_vect_t
    OneD_PeriodicBoundaryCondition::apply(const position_t& ri, const position_t& rj) const {
        static int index = index_(direction_);
        const box_t& box = *box_;

        dist_vect_t rij = ri - rj;
        real_t box_k = box[index];
        real_t dr = ri[index] - rj[index];
        real_t ratio = dr / box_k;
        real_t n = util::nint(ratio);
        rij[index] = dr - n * box_k;

        return rij;
    }

    position_t
    OneD_PeriodicBoundaryCondition::placeInside(const position_t& r_out) const {
        static util::Logger logger("PeriodicBoundaryCondition::placeInside()");
        static int index = index_(direction_);
        const box_t& box = *box_;

        position_t r_in{r_out};
        real_t rk = r_in[index];
        real_t box_k = box[index];
        real_t ratio = rk / box_k;
        real_t n = util::nint(std::floor(ratio));
        r_in[index] -= n * box_k;
        if (r_in[index] < 0 || r_in[index] > box_k) {
            logger.error("r_out:");
            util::log(logger, r_out, util::Logger::LOGERROR);
            logger.error("r_in:");
            util::log(logger, r_in, util::Logger::LOGERROR);
            logger.error("n: " + util::toString(n));
            util::logAndThrow(logger, "Position r_in is -not- inside the box.");
        }
        //assert(r_in[k] >= 0 && r_in[k] <= box_k);

        return r_in;
    }
}

