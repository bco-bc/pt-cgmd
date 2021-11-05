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
 * File:   cg_displacer.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 5, 2019, 12:55 PM
 */

#ifndef CG_DISPLACER_HPP
#define CG_DISPLACER_HPP

#include "displacer.hpp"
#include "sim-data.hpp"
#include "s-types.hpp"

namespace simploce {
    
    /**
     * Displaces beads of a coarse grained simulation model
     */
    struct CoarseGrainedDisplacer : public Displacer {
        
        virtual ~CoarseGrainedDisplacer() {}
        
        /**
         * Displaces beads.
         * @param param Simulation parameters.
         * @param cg Coarse grained particle model.
         * @return Simulation data (e.g. kinetic energy, temperature, etc).
         */
        virtual SimulationData displace(const sim_param_t& param, 
                                        const cg_mod_ptr_t& cg) const = 0;
        
    };
}

#endif /* CG_DISPLACER_HPP */

