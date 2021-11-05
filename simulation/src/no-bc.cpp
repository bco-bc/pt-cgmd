/*
 * File:   no-bc.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:55 PM
 */

#include "simploce/simulation/no-bc.hpp"
#include "simploce/simulation/s-conf.hpp"

namespace simploce {

    dist_vect_t 
    NoBoundaryCondition::apply(const position_t& r1, 
                               const position_t& r2) const 
    {
        return std::move(dist_vect_t{r1 - r2});
    }
    
    position_t 
    NoBoundaryCondition::placeInside(const position_t& r) const
    {
        return r;
    }
    
    std::string 
    NoBoundaryCondition::id() const
    {
        return conf::NOBC;
    }
 
}
