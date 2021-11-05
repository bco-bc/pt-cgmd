/*
 * The MIT License
 *
 * Copyright 2019 juffer.
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
 * File:   langevin-leap-frog.hpp
 * Author: juffer
 *
 * Created on 16 November 2019, 15:23
 */

#ifndef LANGEVIN_LEAP_FROG_HPP
#define LANGEVIN_LEAP_FROG_HPP

#include "cg-displacer.hpp"
#include "simploce/particle/coarse-grained.hpp"

namespace simploce {
    
    /**
     * Displaces particles according to a stochastic impulsive Langevin 
     * Leap-Frog for systems without constrains.
     * @see <a href="https://dx.doi.org/10.1021/ct3000876 | J.">
     *  Goga et al, J. Chem. Theory Comput. 2012, 8, 3637−3649.
     * </a>
     */
    template <typename M>
    class LangevinLeapFrog;
    
    /**
     * Specialization for a coarse grained particle model.
     */
    template <>
    class LangevinLeapFrog<CoarseGrained> : public CoarseGrainedDisplacer {
    public:
        
        LangevinLeapFrog(const cg_interactor_ptr_t& interactor);
        
        /**
         * Returns beads of a coarse grained particle model.
         * @param param Simulation parameters.
         * @param cg Coarse grained particle model.
         * @return Kinetic, potential energy, and temperature.
         */
        SimulationData displace(const sim_param_t& param, 
                                const cg_mod_ptr_t& cg) const override;
        
    private:
        
        cg_interactor_ptr_t interactor_;
    };
}

#endif /* LANGEVIN_LEAP_FROG_HPP */

