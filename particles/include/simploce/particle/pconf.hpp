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
 * File:   conf.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 12, 2019, 4:12 PM
 */

#ifndef PCONF_HPP
#define PCONF_HPP

#include "ptypes.hpp"
#include "simploce/util/uconf.hpp"

namespace simploce {
    namespace conf {
    
        /**
         * Width of a (particle, spec) name output field in an output stream.
         */
        const int NAME_WIDTH = 10;
        
        /**
         * Mass of a proton (u)
         */
        const mass_t MASS_PROTON = 1.007276466879;
        
        /**
         * Charge of a proton (e).
         */
        const charge_t CHARGE_PROTON = 1.0;
        
        /**
         * Discrete changes in particle's protonation state.
         */
        const std::size_t PROTONATABLE = 1;
        
        /**
         * continouos change in particle's protonation state.
         */
        const std::size_t PROTONATING = 2;
        
    }
}

#endif /* CONF_HPP */

