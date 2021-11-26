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
  class PeriodicBoundaryCondition : public boundary_condition {
  public:
      
    /**
     * Constructor.
     * @param box Orthogonal simulation box (container, unit cell). 
     */
    PeriodicBoundaryCondition(box_ptr_t box);

    dist_vect_t apply(const position_t& r1, const position_t& r2) const override;
    
    position_t placeInside(const position_t& r_out) const override;
    
    std::string id() const override;

    void box(const box_ptr_t& box) override;

  private:

    box_ptr_t box_;
    
  };
}


#endif /* PBC_HPP */

