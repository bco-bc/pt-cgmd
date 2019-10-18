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
 * File:   pbc.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:53 PM
 */

#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/stypes.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/util/util.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

namespace simploce {

    PeriodicBoundaryCondition::PeriodicBoundaryCondition(const box_ptr_t& box) :
        BoundaryCondition{}, box_{box}
    {
        if ( !box_ ) {
            throw std::domain_error("PeriodicBoundaryCondition: a box must be provided.");
        }
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
            real_t n = util::nint<int>(dr/boxk);
            std::cout << "boxk, dr, dr/boxk, n " << boxk << ' ' << dr << ' ' 
                      << n << ' ' << dr/boxk << std::endl;
            r12[k] = dr - n * boxk;
        }
        return r12;  
    }
    
    position_t 
    PeriodicBoundaryCondition::placeInside(const position_t& r) const
    {
        const box_t& box = *box_;
        position_t rin{r};
        for (std::size_t k = 0; k != 3; ++k) {
            real_t rk = rin[k];
            real_t boxk = box[k];
            real_t n = std::floor(rk/boxk);
            rin[k] -= n * boxk;
            
            std::cout << "k, r[k], rin[k], n: " << k << ' ' << r[k] << ' ' << rin[k] << ' ' << n << std::endl;
            assert(rin[k] >= 0 && rin[k] <= boxk);
        }
        return rin;        
    }
    
    std::string 
    PeriodicBoundaryCondition::id() const
    {
        return conf::PBC;
    }

}
