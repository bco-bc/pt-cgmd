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
 * File:   langevin-velocity-verlet.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 9:47 AM
 */

#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/sproperties.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/util/mu-units.hpp"
#include "simploce/util/util.hpp"
#include <random>
#include <cmath>
#include <array>

namespace simploce {
    
    // Helpers.
    static std::vector<real_t> frictionCoefficients_{};
    static std::vector<real_t> B_{};
    static std::vector<real_t> A1_{};
    static std::vector<real_t> A2_{};
    static std::vector<real_t> strengths_{};
    
    static std::vector<force_t> fis_{};      // Forces at time t(n).
    static std::vector<position_t> ris_{};   // Positions at time t(n).
    static std::vector<momentum_t> pis_{};   // Momenta at time t(n).
    
    // Random vector, each element is a array of size 3.
    static std::random_device rd_;
    static std::mt19937 gen_(rd_());
    static std::normal_distribution<real_t> dis_;
    static std::vector<std::array<real_t, 3>> W_{};
    
    /**
     * Calculation of helper values.
     * @param T Particle type.
     * @param dt Time step.
     * @param particles Particles.
     */
    template <typename T>
    void setupHelpers_(const stime_t& dt,
                       const temperature_t& temperature,
                       real_t gamma,
                       const std::vector<std::shared_ptr<T>>& particles)
    {
        real_t kT = MUUnits<real_t>::KB * temperature();
        std::size_t nparticles = particles.size();
        
        frictionCoefficients_ = std::vector<real_t>(nparticles, 0.0);
        B_ = std::vector<real_t>(nparticles, 0.0);
        A1_ = std::vector<real_t>(nparticles, 0.0);
        A2_ = std::vector<real_t>(nparticles, 0.0);
        strengths_ = std::vector<real_t>(nparticles, 0.0);
        
        fis_ = std::vector<force_t>(nparticles, force_t{});
        ris_ = std::vector<position_t>(nparticles, position_t{});
        pis_ = std::vector<momentum_t>(nparticles, momentum_t{});
        
        W_ = std::vector<std::array<real_t, 3>>(nparticles, std::array<real_t, 3>{0.0, 0.0, 0.0});
        auto value = util::seedValue<std::size_t>();
        gen_.seed(value);
    
        for (std::size_t index = 0; index != nparticles; ++index) {
            auto& particle = *particles[index];
        
            mass_t mass = particle.mass();                        // In u.
            mass_t mass2 = mass * 2.0;
            real_t fc =  mass() * gamma;                          // Friction coefficient in u/ps.
            frictionCoefficients_[index] = fc;                    // Save friction coefficient for 
                                                                  // momentum update.

            // All kinds of constant factors for each particle.
            B_[index] = 1.0 / ( 1.0 + fc * dt() / mass2() );      // No units.
            real_t a1 = dt() / mass2();                           // ps/u
            A1_[index] = a1;
            real_t a2 = a1 * dt();                                // ps^2/u
            A2_[index] = a2;
            real_t strength = std::sqrt( dt() * 2.0 * fc * kT );  // In (u nm) / ps.
            strengths_[index] = strength;

            // Validate.
            real_t f1 = fc * a1;
            if ( f1 > 1) {
                std::clog << "WARNING: " << std::endl;
                std::clog << "Langevin Velocity Verlet: " << std::endl;
                std::clog << "f1 = " << f1 << ": f1 > 1 may represent an unphysical regime." << std::endl;
            }
        }                
    }

    /**
     * Displace particle position.
     * @param T Particle type.
     * @param dt Time step.
     * @param particles Particles.
     */
    template <typename T>
    void displacePosition_(const stime_t& dt,
                           const std::vector<std::shared_ptr<T>>& particles)
    {        
        for (std::size_t index = 0; index != particles.size(); ++index) {
            auto& particle = *particles[index];
    
            // Some required values.
            mass_t mass = particle.mass();

            // Update position, not momentum.
            static std::array<real_t, 3> w;            // Random vector at t(n+1).
            w[0] = dis_(gen_);
            w[1] = dis_(gen_);
            w[2] = dis_(gen_);
            W_[index] = w;                             // Save for updating velocities.
        
            force_t fi = particle.force();             // Force (kJ/(mol nm) = (u nm)/(ps^2))
            fis_[index] = fi;                          // at time t(n).
      
            momentum_t pi = particle.momentum();       // Momentum (u nm/ps) at time t(n).
            pis_[index] = pi;                          // Save for next step.
            position_t ri = particle.position();       // Position at time t(n).
            ris_[index] = ri;                          // Save for momentum update.
            static position_t rf;
            real_t b = B_[index];                      // No units.
            real_t a1 = A1_[index];                    // ps/u
            real_t a2 = A2_[index];                    // ps^2/u
            real_t strength = strengths_[index];
            for (std::size_t k = 0; k != 3; ++k) { 
                real_t vi = pi[k] / mass();            // Velocity (nm/ps) at time t(n)
                rf[k] = ri[k] +
                    b * dt() * vi +
                    b * a2 * fi[k] +
                    b * a1 * strength * w[k];          // Position at time t(n+1).
            }
        
            // Update position.
            particle.position(rf);        
        }
    }
    
    template <typename T>
    SimulationData displaceMomentum_(const std::vector<std::shared_ptr<T>>& particles)
    {
        SimulationData data;
        
        // Kinetic energy at t(n+1).
        data.ekin = 0.0;
        
        // Displace particles: momenta.        
        for (std::size_t index = 0; index != particles.size(); ++index) {
            auto& particle = *particles[index];
            
            mass_t mass = particle.mass();             // In u.
      
            std::array<real_t, 3> w = W_[index];       // Random vector at t(n+1).
            force_t fi = fis_[index];                  // Force (kJ/(mol nm) = (u nm)/(ps^2)) 
                                                       // at time t(n).
            position_t ri = ris_[index];               // Position at time t(n).
      
            force_t ff = particle.force();             // Force (kJ/(mol nm) = (u nm)/(ps^2))
                                                       // at time t(n+1).
            momentum_t pi = particle.momentum();       // Momentum (u nm/ps) at time t(n).
            position_t rf = particle.position();       // Position at time t(n+1).
      
            real_t fc = frictionCoefficients_[index];
            real_t a1 = A1_[index];
            real_t strength = strengths_[index];
      
            static velocity_t vf{};
            static momentum_t pf{};
            for (std::size_t k = 0; k != 3; ++k) {
                vf[k] = pi[k] / mass() +
                    a1 * ( fi[k] + ff[k] ) -
                    fc * ( rf[k] - ri[k] ) / mass() +
                    strength * w[k] / mass();          // Velocity at time t(n+1).
                pf[k] = mass() * vf[k];                // Momentum (u nm/ps) at time t(n+1). 
            }
      
            // Update momentum.
            particle.momentum(pf);
      
            // Kinetic energy at time t(n+1)..
            data.ekin += 0.5 * mass() * inner<real_t>(vf, vf);
        }
        
        // Instantaneous temperature at t(n+1).
        data.temperature = temperature<T>(particles, data.ekin);

        // Done.
        return data;        
        
    }
    
    LangevinVelocityVerlet<Atomistic>::LangevinVelocityVerlet(const at_interactor_ptr_t& interactor) :
        AtomisticDisplacer{}, interactor_{interactor}
    {       
    }
        
    SimulationData 
    LangevinVelocityVerlet<Atomistic>::displace(const sim_param_t& param, 
                                                const at_ptr_t& at) const
    {
        static bool setup = false;

        static stime_t dt{0.0};
        static temperature_t temperature{298.15};
        static real_t gamma{0.5};
        static std::size_t counter = 0;

        if ( !setup ) {            
            dt = param.get<real_t>("timestep");
            temperature = param.get<real_t>("temperature", 298.15);
            gamma = param.get<real_t>("gamma", 0.5);

            at->doWithAll<void>([] (const std::vector<atom_ptr_t>& atoms) {
                setupHelpers_<Atom>(dt, temperature, gamma, atoms);
            });
            
            interactor_->interact(param, at);  // Initial forces.
            
            setup = true;
        }
        
        counter += 1;
                
        // Displace atom positions.
        at->doWithAll<void>([] (const std::vector<atom_ptr_t>& atoms) {
            displacePosition_<Atom>(dt, atoms);
        });
        
        // Compute forces and potential energy at t(n+1) using positions at t(n+1).
        energy_t epot = interactor_->interact(param, at);
        
        SimulationData data = at->doWithAll<SimulationData>([] (const std::vector<atom_ptr_t>& atoms) {
            return displaceMomentum_<Atom>(atoms);
        });
        data.epot = epot;
        
        data.t = counter * dt;
        
        return data;
    }    
    
    std::string LangevinVelocityVerlet<Atomistic>::id() const
    {
        return conf::LANGEVIN_VELOCITY_VERLET;
    }

    LangevinVelocityVerlet<CoarseGrained>::LangevinVelocityVerlet(const cg_interactor_ptr_t& interactor) :
        CoarseGrainedDisplacer{}, interactor_{interactor}     
    {       
    }

    SimulationData 
    LangevinVelocityVerlet<CoarseGrained>::displace(const sim_param_t& param, 
                                                    const cg_ptr_t& cg) const
    {
        static bool setup = false;

        static stime_t dt{0.0};
        static temperature_t temperature{0.0};
        static real_t gamma{0.0};
        static std::size_t counter = 0.0;
        
        counter += 1;
        
        if ( !setup ) {
            dt = param.get<real_t>("timestep");
            temperature = param.get<real_t>("temperature");
            gamma = param.get<real_t>("gamma");
            
            cg->doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
                setupHelpers_<Bead>(dt, temperature, gamma, beads);
            });
            
            interactor_->interact(param, cg); // Initial forces.
            
            setup = true;
        }
        
        // Displace atom positions.
        cg->doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
            displacePosition_<Bead>(dt, beads);
        });
        
        // Compute forces and potential energy at t(n+1) using positions at t(n+1).
        energy_t epot = interactor_->interact(param, cg);
        
        SimulationData data = cg->doWithAll<SimulationData>([] (const std::vector<bead_ptr_t>& beads) {
            return displaceMomentum_<Bead>(beads);
        });
        data.epot = epot;
        
        data.t = counter * dt;

        return data;
    }    

    std::string LangevinVelocityVerlet<CoarseGrained>::id() const
    {
        return conf::LANGEVIN_VELOCITY_VERLET;
    }

}

