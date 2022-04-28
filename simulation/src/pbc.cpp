/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:53 PM
 */

#include "simploce/simulation/pbc.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

namespace simploce {

    PeriodicBoundaryCondition::PeriodicBoundaryCondition(box_ptr_t box) :
            boundary_condition{}, box_{std::move(box)} {
    }

    dist_vect_t 
    PeriodicBoundaryCondition::apply(const position_t& ri,
                                     const position_t& rj) const
    {        
        const box_t& box = *box_;
    
        dist_vect_t rij{};
        for (std::size_t k = 0; k != 3; ++k) {
            real_t box_k = box[k];
            real_t dr = ri[k] - rj[k];
            real_t ratio = dr / box_k;
            real_t n = util::nint(ratio);
            rij[k] = dr - n * box_k;
        }
        return rij;
    }
    
    position_t 
    PeriodicBoundaryCondition::placeInside(const position_t& r_out) const
    {
        static util::Logger logger("PeriodicBoundaryCondition::placeInside()");
        const box_t& box = *box_;
        
        position_t r_in{r_out};
        for (std::size_t k = 0; k != 3; ++k) {
            real_t rk = r_in[k];
            real_t box_k = box[k];
            real_t ratio = rk / box_k;
            real_t n = util::nint(std::floor(ratio));
            r_in[k] -= n * box_k;
            if (r_in[k] < 0 || r_in[k] > box_k) {
                logger.error("r_out:");
                util::log(logger, r_out, util::Logger::LOGERROR);
                logger.error("r_in:");
                util::log(logger, r_in, util::Logger::LOGERROR);
                logger.error("n: " + util::toString(n));
                util::logAndThrow(logger, "Position r_in is -not- inside the box.");
            }
            //assert(r_in[k] >= 0 && r_in[k] <= box_k);
        }
        return r_in;
    }

}
