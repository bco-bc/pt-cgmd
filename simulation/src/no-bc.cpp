/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:55 PM
 */

#include "simploce/simulation/no-bc.hpp"
#include "simploce/simulation/s-conf.hpp"

namespace simploce {

    dist_vect_t 
    NoBoundaryCondition::apply(const position_t& ri,
                               const position_t& rj) const
    {
        return std::move(dist_vect_t{ri - rj});
    }
    
    position_t 
    NoBoundaryCondition::placeInside(const position_t& r_out) const
    {
        return r_out;
    }
    
    std::string 
    NoBoundaryCondition::id() const
    {
        return conf::NO_BC;
    }
 
}
