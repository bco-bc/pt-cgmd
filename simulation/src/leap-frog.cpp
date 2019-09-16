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
 * File:   leap-frog.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 5:12 PM
 */

#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/sproperties.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/bead.hpp"
#include <vector>

namespace simploce {
    
    /*
     * @param T particle pointer type.
     * @param dt Time step. 
     * @param particles Particles.
     * @return Kinetic energy, temperature.
     */
    template <typename T>
    static SimulationData displace_(const stime_t dt, 
                                    const std::vector<std::shared_ptr<T>>& particles)
    {
        
        // Assume current step n-1/2 at time t(n-1/2).
            
        // Compute linear momentum and position, plus kinetic energy.
        SimulationData data;
        for (auto ptr : particles) {
            T &particle = *ptr;
            mass_t mass = particle.mass();             // In u.
            const force_t &f = particle.force();       // Force (kJ/(mol nm) at time t(n-1/2).
            momentum_t p = particle.momentum();        // Momentum (u nm/ps) at time t(n-1/2).
        
            static velocity_t vi{}, vf{};
            position_t r = particle.position();        // Position at time t(n).
      
            for (std::size_t k = 0; k != 3; ++k) {
                vi[k] = p[k] / mass();                     // Velocity (nm/ps) at time t(n-1/2)
                p[k] += dt() * f[k];                      // Momentum at time t(n+1/2)
        
                vf[k] = p[k] / mass();                     // Velocity at time t(n+1/2)
                r[k] += dt() * vf[k];                     // Position at time t(n+1).
            }

            // Save new position and linear momentum.
            particle.position(r);
            particle.momentum(p);

            // Kinetic energy
            velocity_t va = 0.5 * (vi + vf);           // Average velocity at time t(n).
            data.ekin += 0.5 * mass() * inner<real_t>(va, va);
        }
        
        // Temperature at t(n).
        data.temperature = temperature<T>(particles, data.ekin);
        
        return data;        
    }
    
    LeapFrog<Atomistic>::LeapFrog(const at_interactor_ptr_t& interactor) :
        interactor_{interactor}
    {        
    }
    
    SimulationData LeapFrog<Atomistic>::displace(const sim_param_t& param, 
                                                 const at_ptr_t& at) const
    {
        static bool setup = false;
        static stime_t dt{0.0};
        
        if ( !setup ) {
            dt = param.get<real_t>("timestep");        
        }
        
        interactor_->interact(param, at);
        return at->doWithAll<SimulationData>([] (const std::vector<atom_ptr_t>& atoms) {
            return displace_<Atom>(dt, atoms);
        });        
    }
    
    std::string LeapFrog<Atomistic>::id() const
    {
        return conf::LEAP_FROG;
    }
    
    LeapFrog<CoarseGrained>::LeapFrog(const cg_interactor_ptr_t& interactor) : 
        interactor_{interactor}
    {        
    }
    
    SimulationData LeapFrog<CoarseGrained>::displace(const sim_param_t& param, 
                                                     const cg_ptr_t& cg) const
    {
        static bool setup = false;
        static stime_t dt{0.0};
        
        if ( !setup ) {
            dt = param.get<real_t>("timestep");        
        }
        
        interactor_->interact(param, cg);
        return cg->doWithAll<SimulationData>([] (const std::vector<bead_ptr_t>& beads) {
            return displace_<Bead>(dt, beads);
        });
    }
    
    std::string LeapFrog<CoarseGrained>::id() const
    {
        return conf::LEAP_FROG;
    }
}

