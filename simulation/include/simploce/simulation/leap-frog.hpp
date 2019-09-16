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
 * File:   leap-frog.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:36 PM
 */

#ifndef LEAP_FROG_HPP
#define LEAP_FROG_HPP

#include "cg-displacer.hpp"
#include "at-displacer.hpp"
#include "sim-data.hpp"
#include "stypes.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"

namespace simploce {        
    
    /**
     * Standard leap frog algorithm for MD simulations.
     * @param M Particle model type.
     */
    template <typename M>
    class LeapFrog;
    
    /**
     * Specialization for an atomistic particle model.
     */
    template <>
    class LeapFrog<Atomistic> : public AtomisticDisplacer{
    public:    
        
        LeapFrog(const at_interactor_ptr_t& interactor);
                
        /**
         * Displaces atoms of an atomistic model.
         * @param at Atomistic model.
         * @return Kinetic energy, temperature.
         */
        SimulationData displace(const sim_param_t& param, 
                                const at_ptr_t& at) const override;
        
        std::string id() const override;
        
    private:
        
       at_interactor_ptr_t interactor_;
        
    };
    
    /**
     * Specialization for a coarse grained particle model.
     */
    template <>
    class LeapFrog<CoarseGrained> : public CoarseGrainedDisplacer {
    public:
        
        LeapFrog(const cg_interactor_ptr_t& interactor);
        
        SimulationData displace(const sim_param_t& param, 
                                const cg_ptr_t& cg) const override;
        
        std::string id() const override;
        
    private:
    
        cg_interactor_ptr_t interactor_;
    };
    
        
}

#endif /* LEAP_FROG_HPP */

