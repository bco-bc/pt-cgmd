/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:53 PM
 */

#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/util/util.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

namespace simploce {

    PeriodicBoundaryCondition::PeriodicBoundaryCondition(const box_ptr_t& box) :
            boundary_condition{}, box_{factory::box(box->edgeLength())} {
    }

    dist_vect_t 
    PeriodicBoundaryCondition::apply(const position_t& r1, 
                                     const position_t& r2) const 
    {        
        const box_t& box = *box_;
    
        dist_vect_t r12{};
        for (std::size_t k = 0; k != 3; ++k) {
            real_t boxk = box[k];
            real_t dr = r1[k] - r2[k];
            real_t ratio = dr/boxk;
            real_t n = util::nint(ratio);
            r12[k] = dr - n * boxk;
        }
        return r12;  
    }
    
    position_t 
    PeriodicBoundaryCondition::placeInside(const position_t& r_out) const
    {
        const box_t& box = *box_;
        
        position_t rin{r_out};
        for (std::size_t k = 0; k != 3; ++k) {
            real_t rk = rin[k];
            real_t boxk = box[k];
            real_t n = util::nint(std::floor(rk/boxk));
            rin[k] -= n * boxk;
            assert(rin[k] >= 0 && rin[k] <= boxk);
        }
        return rin;        
    }
    
    std::string 
    PeriodicBoundaryCondition::id() const
    {
        return conf::PBC;
    }

    void
    PeriodicBoundaryCondition::box(const box_ptr_t& box) {
        box_ = factory::box(box->edgeLength());
    }

}
