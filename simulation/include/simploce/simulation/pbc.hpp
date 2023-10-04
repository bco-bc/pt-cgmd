/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:50 PM
 */

#ifndef PBC_HPP
#define PBC_HPP

#include "bc.hpp"

namespace simploce {

  /**
   * Periodic boundary condition with the nearest image approximation.
   * @see <a href="https://en.wikipedia.org/wiki/Periodic_boundary_conditions">Wikipedia</a>
   */
  class PBC : public boundary_condition {
  public:
      
    /**
     * Constructor.
     * @param box Orthogonal simulation box (container, unit findCell).
     */
    explicit PBC(box_ptr_t box);

    dist_vect_t apply(const position_t& ri, const position_t& rj) const override;
    
    position_t placeInside(const position_t& r_out) const override;

    /**
     * No effect.
     */
    velocity_t apply(const simploce::velocity_t &v,
                     const position_t& r) const override;

    /**
     * No effect.
     */
    void applyToVelocities(const pg_ptr_t& particleGroup) const override;
    
  private:

    box_ptr_t box_;
    
  };
}


#endif /* PBC_HPP */

