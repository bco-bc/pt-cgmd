/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:49 PM
 */

#ifndef BC_HPP
#define BC_HPP

#include "s-types.hpp"

namespace simploce {

  /**
   * As applicable to a given molecular simulation for calculating interactions.
   */
  struct boundary_condition {

    virtual ~boundary_condition() = default;

    /**
     * Applies the boundary condition.
     * @param r1 Position.
     * @param r2 Position.
     * @return Difference r1-r2 under the given boundary condition.
     */
    virtual dist_vect_t apply(const position_t& r1, const position_t& r2) const = 0;
    
    /**
     * Moves given position to the inside of the simulation box.
     * @param r_out Position possibly outside simulation box.
     * @return Position Inside simulation box.
     */
    virtual position_t placeInside(const position_t& r_out) const = 0;
    
    /**
     * Returns an identifier.
     * @return Identifier.
     */
    virtual std::string id() const = 0;

    /**
     * Sets or updates the box.
     * @param box.
     */
    virtual void box(const box_ptr_t& box) {};

  };
  
}


#endif /* BC_HPP */

