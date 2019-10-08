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
 * File:   forcefield.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 13, 2019, 3:22 PM
 */

#ifndef FORCEFIELD_HPP
#define FORCEFIELD_HPP

#include "stypes.hpp"
#include "simploce/util/map2.hpp"
#include <string>
#include <utility>

namespace simploce {
    
    /**
     * Force field.
     */
    struct ForceField {
        
        /**
         * Type for holding pairs (C12, C6) of LJ interaction parameters 
         * for pairs of particle specs, each identified by a specification name. 
         */
        using lj_params_t = MatrixMap<std::string, std::pair<real_t, real_t>>;
        
        /**
         * Type for holding parameters for electrostatic interactions.
         */
        using el_params_t = std::map<std::string, real_t>;
        
        virtual ~ForceField() {}
        
        /**
         * Returns identifying name.
         * @return Identifying name.
         */
        virtual std::string id() const = 0;
        
        /**
         * Returns Lennard-Jones and electrostatic interaction parameters.
         * @return Parameters,
         */
        virtual std::pair<lj_params_t, el_params_t> parameters() const = 0;
    };
}

#endif /* FORCEFIELD_HPP */

