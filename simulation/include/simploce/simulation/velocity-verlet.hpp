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
 * File:   velocity-verlet.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:38 PM
 */

#ifndef VELOCITY_VERLET_HPP
#define VELOCITY_VERLET_HPP

#include "cg-displacer.hpp"
#include "at-displacer.hpp"
#include "sim-data.hpp"
#include "stypes.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"

namespace simploce {
    
    /**
     * Velocity Verlet algorithm for MD simulations.
     * @param M Particle model type.
     */
    template <typename M>
    class VelocityVerlet;
    
    /**
     * Specialization for atomistic model.
     */
    template <>
    class VelocityVerlet<Atomistic> : public AtomisticDisplacer {
    public:
        
        VelocityVerlet(const at_interactor_ptr_t& interactor);
        
        /**
         * Displaces beads of an atomistic model.
         * @param at Atomistic model.
         * @return kinetic, potential energy, and temperature.
         */
        SimulationData displace(const sim_param_t& param, 
                                const at_ptr_t& at) const override;
        
        std::string id() const override;
        
    private:
        
        at_interactor_ptr_t interactor_;
    };
    
    /**
     * Specialization for coarse grained model.
     */
    template <>
    class VelocityVerlet<CoarseGrained> : public CoarseGrainedDisplacer {
    public:
        
        VelocityVerlet(const cg_interactor_ptr_t& interactor);
        
        /**
         * Displaces beads of an atomistic model.
         * @param at Atomistic model.
         * @return kinetic, potential energy, and temperature.
         */
        SimulationData displace(const sim_param_t& param, 
                                const cg_ptr_t& cg) const override;
        
        std::string id() const override;
        
    private:
        
        cg_interactor_ptr_t interactor_;
            
    };
    
}

#endif /* VELOCITY_VERLET_HPP */

