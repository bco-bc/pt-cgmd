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
 * File:   sconf.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 11:37 AM
 */

#ifndef SCONF_HPP
#define SCONF_HPP

#include "simploce/particle/pconf.hpp"
#include "stypes.hpp"

namespace simploce {
    namespace conf {
        
        const std::string LEAP_FROG = "lf";
        const std::string LANGEVIN_VELOCITY_VERLET = "lvv";
        const std::string PT_LANGEVIN_VELOCITY_VERLET = "pt-lvv";
        const std::string VELOCITY_VERLET = "vv";
        
        const std::string NOBC = "no-bc";
        const std::string PBC = "pbc";
        
        const std::string POLARIZABLE_WATER = "pol-water";
        const std::string ACID_BASE_SOLUTION = "acid-base-solution";
        const std::string ELECTROLYTE = "electrolyte";
        
        // Default cutoff distance for non bonded interactions.
        static length_t RCUTOFF_DISTANCE_{2.5};  // nm.
        
        // Minimum number of particles.
        const std::size_t MIN_NUMBER_OF_PARTICLES = 1000;    
        
        /**
         * Default cutoff distance beyond which protonatables cannot transfer protons.
         */
        static length_t RCUTOFF_DISTANCE_PT{0.4};
    }
}

#endif /* SCONF_HPP */

