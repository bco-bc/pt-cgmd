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
 * File:   properties.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 3:24 PM
 */

#ifndef SPROPERTIES_HPP
#define SPROPERTIES_HPP

#include "simploce/util/mu-units.hpp"

namespace simploce {
    
    /**
     * Calculates instantaneous temperature for a collection of particles.
     * @param particles Particles.
     * @param ekin Kinetic energy.
     * @return Instantaneous temperature.
     */
    template <typename T>
    temperature_t temperature(const std::vector<std::shared_ptr<T>>& particles, 
                              const energy_t& ekin)
    {
        std::size_t nparticles = particles.size();
        //std::size_t ndof = ( 3 * nparticles - 3 );  // -3 to remove rigid body translation
        std::size_t ndof = 3 * nparticles;
        if ( ndof > 3 ) {
            return 2.0 * ekin() / ( real_t(ndof) * MUUnits<real_t>::KB );  // In K.
        } else {
            // No point calculating temperature for a low number of degrees of freedom.
            return 0.0;
        }        
    }
    
}

#endif /* PROPERTIES_HPP */

