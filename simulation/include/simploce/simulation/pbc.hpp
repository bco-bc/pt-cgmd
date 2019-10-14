/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   pbc.hpp
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
  class PeriodicBoundaryCondition : public BoundaryCondition {
  public:
      
    /**
     * Constructor.
     * @param box Orthogonal simulation box (container, unit cell). 
     */
    PeriodicBoundaryCondition(const box_ptr_t& box);

    dist_vect_t apply(const position_t& r1, const position_t& r2) const override;
    
    position_t placeInside(const position_t& r) const override;
    
    std::string id() const override;

  private:

    box_ptr_t box_;
    
  };
}


#endif /* PBC_HPP */

