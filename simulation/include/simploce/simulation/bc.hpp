/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:49 PM
 */

#ifndef BC_HPP
#define BC_HPP

#include "simploce/types/s-types.hpp"

namespace simploce {

  /**
   * As applicable to a given molecular simulation for calculating interactions.
   */
  struct boundary_condition {

    virtual ~boundary_condition() = default;

    /**
     * Applies the boundary condition.
     * @param ri Position.
     * @param rj Position.
     * @return Difference ri-rj under the given boundary condition.
     */
    virtual dist_vect_t apply(const position_t& ri, const position_t& rj) const = 0;
    
    /**
     * Moves given position inside simulation box.
     * @param r_out Position possibly outside simulation box.
     * @return Position Inside simulation box.
     */
    virtual position_t placeInside(const position_t& r_out) const = 0;

  };
  
}


#endif /* BC_HPP */

