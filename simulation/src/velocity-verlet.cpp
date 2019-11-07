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
 * File:   velocity-verlet.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 5:39 PM
 */

#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/sim-util.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/bead.hpp"
#include <utility>

namespace simploce {
    
    /*
     * Displaces particle positions.
     * @param T Particle type.
     * @param dt Time step.
     * @param particles Particles.
     * @return Forces at time t(n).
     */
    template <typename T>
    std::vector<force_t> 
    displacePosition_(const stime_t& dt,
                      const std::vector<std::shared_ptr<T>>& particles)
    {
        // Initial forces (forces at time t(n).
        std::vector<force_t> fis(particles.size(), force_t{});
        
        // Displace particles: Positions.
        for (std::size_t index = 0; index != particles.size(); ++index) {
            T& particle = *particles[index];

            mass_t mass = particle.mass();
            real_t a1 = dt() / ( 2.0 * mass() );
            real_t a2 = dt() * a1;

            force_t fi = particle.force();                 // Force (kJ/(mol nm) = (u nm)/(ps^2)) 
            fis[index] = fi;                               // at time t(n).
      
            position_t ri = particle.position();           // Position at time t(n).
            velocity_t vi = particle.velocity();           // Velocity at time t(n).
            static position_t rf{};
            for ( std::size_t k = 0; k != 3; ++k) {
                rf[k] = ri[k] + dt() * vi[k] + a2 * fi[k]; // Position at time t(n+1).
            }
        
            // Save new position.       
            particle.position(rf);
        }
        return std::move(fis);
    }
    
    /*
     * Displaces particle velocities.
     * @param T Particle type.
     * @param dt Time step.
     * @param fis Forces at t(n)
     * @param particles Particles.
     * @return Kinetic, potential energy, and temperature.
     */
    template <typename T>
    SimulationData 
    displaceMomentum_(const stime_t& dt,
                      const std::vector<force_t>& fis,
                      const std::vector<std::shared_ptr<T>>& particles)
    {
        SimulationData data;
        
        // Kinetic energy at t(n+1).
        data.ekin = 0.0;
        
        // Displace particles: Momenta/velocities.
        std::size_t nparticles = particles.size();
        for (std::size_t index = 0; index != nparticles; ++index) {
            T& particle = *particles[index];
            mass_t mass = particle.mass();
            real_t a1 = dt() / ( 2.0 * mass() );

            force_t fi = fis[index];                       // Force (kJ/(mol nm) = (u nm)/(ps^2)) 
                                                           // at time t(n).
            force_t ff = particle.force();                 // Force (kJ/(mol nm) = (u nm)/(ps^2))
                                                           // at time t(n+1).

            velocity_t vi = particle.velocity();           // velocity (nm/ps) at time t(n).
            static velocity_t vf{};
            for (std::size_t k = 0; k != 3; ++k) {
                vf[k] = vi[k] + a1 * ( fi[k] + ff[k] );    // Velocity at time t(n+1).
            }
      
            // Save velocity
            particle.velocity(vf);
      
            // Kinetic energy at t(n+1).
            data.ekin += 0.5 * mass() * inner<real_t>(vf, vf);
        }
    
        // Instantaneous temperature at t(n+1).
        data.temperature = util::temperature<T>(particles, data.ekin);
        
        return data;        
    }
       
    VelocityVerlet<Atomistic>::VelocityVerlet(const at_interactor_ptr_t& interactor) :
        interactor_{interactor}
    {       
    }
        
    SimulationData 
    VelocityVerlet<Atomistic>::displace(const sim_param_t& param, 
                                        const at_ptr_t& at) const
    {        
        static std::size_t counter = 0;
        static bool firstTime = true;
        
        static stime_t dt{0.0};        
        static std::vector<force_t> fis{};
        
        counter += 1;
        
        if ( firstTime) {
            dt = param.get<real_t>("timestep");
            interactor_->interact(param, at); // Initial forces.
            firstTime = false;
        }
        
        // Displace atom positions.
        at->doWithAll<void>([] (const std::vector<atom_ptr_t>& atoms) {
            fis = displacePosition_<Atom>(dt, atoms);
        });
        
        // Compute forces and potential energy at t(n+1) using positions at t(n+1).
        auto result = interactor_->interact(param, at);
        
        // Displace atom momenta.
        SimulationData data = at->doWithAll<SimulationData>([] (const std::vector<atom_ptr_t>& atoms) {
            return displaceMomentum_<Atom>(dt, fis, atoms);
        });
        
        // Save simulation data.
        data.bepot = result.first;
        data.nbepot = result.second;       
        data.t = counter * dt;
        
        return data;
    }
    
    std::string 
    VelocityVerlet<Atomistic>::id() const
    {
        return conf::VELOCITY_VERLET;
    }
        
    VelocityVerlet<CoarseGrained>::VelocityVerlet(const cg_interactor_ptr_t& interactor) :
        interactor_{interactor}
    {       
    }
        
    SimulationData 
    VelocityVerlet<CoarseGrained>::displace(const sim_param_t& param, 
                                            const cg_ptr_t& cg) const
    {        
        static std::size_t counter = 0;
        static bool firstTime = true;
        
        static stime_t dt{0.0};      
        static std::vector<force_t> fis{};
        
        counter += 1;
        if ( firstTime) {
            dt = param.get<real_t>("timestep");
            interactor_->interact(param, cg);
            firstTime = false;
        }
        
        // Displace atom positions.
        cg->doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
            fis = displacePosition_<Bead>(dt, beads);
        });
        
        // Compute forces and potential energy at t(n+1) using positions at t(n+1).
        auto result = interactor_->interact(param, cg);
        
        // Displace atom momenta.
        SimulationData data = cg->doWithAll<SimulationData>([] (const std::vector<bead_ptr_t>& beads) {
            return displaceMomentum_<Bead>(dt, fis, beads);
        });
        
        // Save simulation data.
        data.bepot = result.first;
        data.nbepot = result.second;        
        data.t = counter * dt;
        
        return data;
    }
        
    std::string 
    VelocityVerlet<CoarseGrained>::id() const
    {
        return conf::VELOCITY_VERLET;
    }
        
    
}
