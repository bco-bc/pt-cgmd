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
 * File:   langevin-velocity-verlet.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 4:04 PM
 */

#ifndef LANGEVIN_VELOCITY_VERLET_HPP
#define LANGEVIN_VELOCITY_VERLET_HPP

#include "cg-displacer.hpp"
#include "at-displacer.hpp"
#include "sim-data.hpp"
#include "stypes.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"

namespace simploce {
    
    
    /**
     * Displaces particles according a stochastic Verlet-type algorithm applicable to 
     * an Langevin equation. Provides for a canonical ensemble (NVT constant) simulation. 
     * Requires force field, temperature, time step and damping rate.
     * @see <a href="http://dx.doi.org/10.1080/00268976.2012.760055">
     *   Grønbech-Jensen and Oded Farago, Molec Phys,111, 983-991, 2013
     * </a>
     * @param M Particle model type.
    */
    template <typename M>
    class LangevinVelocityVerlet;
    
    /**
     * Specialization for atomistic particle model.
     */
    template <>
    class LangevinVelocityVerlet<Atomistic> : public AtomisticDisplacer {
    public:
 
        LangevinVelocityVerlet(const at_interactor_ptr_t& interactor);
        
        
        /**
         * Displaces atoms of an atomistic model.
         * @param at Atomistic particle model.
         * @return kinetic, potential energy, and temperature.
         */
        SimulationData displace(const sim_param_t& param, 
                                const at_ptr_t& at) const override;
        
        std::string id() const override;
                
    private:
        
        at_interactor_ptr_t interactor_;
        
    };
    
    
    /**
     * Specialization for coarse grained particle model.
     */
    template <>
    class LangevinVelocityVerlet<CoarseGrained> : public CoarseGrainedDisplacer {
    public:
        
        LangevinVelocityVerlet(const cg_interactor_ptr_t& interactor);
        
        
        /**
         * Displaces atoms of an atomistic model.
         * @param cg Coarse grained particle model.
         * @return kinetic, potential energy, and temperature.
         */
        SimulationData displace(const sim_param_t& param, 
                                const cg_ptr_t& cg) const override;
        
        std::string id() const override;
                
    private:
        
        cg_interactor_ptr_t interactor_;
        
    };
    
}

#endif /* LANGEVIN_VELOCITY_VERLET_HPP */

