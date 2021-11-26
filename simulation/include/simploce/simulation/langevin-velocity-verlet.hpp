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

#include "displacer.hpp"

namespace simploce {
    
    
    /**
     * Displaces particles according a stochastic Verlet-typeName algorithm applicable to
     * an Langevin equation. Provides for a canonical ensemble (NVT constant) simulation. 
     * Requires force field, temperature, time step and damping rate.
     * @see <a href="http://dx.doi.org/10.1080/00268976.2012.760055">
     *   Grønbech-Jensen and Oded Farago, Molec Phys,111, 983-991, 2013
     * </a>
     * @param M Particle model typeName.
    */
    class LangevinVelocityVerlet : public Displacer {
    public:
 
        LangevinVelocityVerlet(sim_param_ptr_t simulationParameters,
                               interactor_ptr_t interactor);

        SimulationData displace(const p_system_ptr_t& particles) const override;

    private:

        sim_param_ptr_t simulationParameters_;
        interactor_ptr_t interactor_;
        
    };

}

#endif /* LANGEVIN_VELOCITY_VERLET_HPP */

