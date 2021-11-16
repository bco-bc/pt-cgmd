/*
 * File:   leap-frog.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 5:12 PM
 */

#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/s-conf.hpp"
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
    static SimulationData 
    displace_(const stime_t dt, 
              const std::vector<std::shared_ptr<T>>& particles)
    {
        static std::size_t counter = 0;
        
        // Assume current step n-1/2 at time t(n-1/2).

        counter += 1;
        
        // Compute linear momentum and position, plus kinetic energy.
        SimulationData data;
        for (auto ptr : particles) {
            T &particle = *ptr;
            mass_t mass = particle.mass();             // In u.
            const force_t &f = particle.force();       // Force (kJ/(mol nm) at time t(n-1/2).
            velocity_t vi = particle.velocity();       // velocity (nm/ps) at time t(n-1/2).
            position_t r = particle.position();        // Position at time t(n).
      
            static velocity_t vf{};
            for (std::size_t k = 0; k != 3; ++k) {
                vf[k] = vi[k] + dt() * f[k] / mass();  // Velocity at time t(n+1/2)
                r[k] += dt() * vf[k];                  // Position at time t(n+1).
            }

            // Save new position and velocity.
            particle.position(r);
            particle.velocity(vf);

            // Kinetic energy
            velocity_t va = 0.5 * (vi + vf);           // Average velocity at time t(n).
            data.ekin += 0.5 * mass() * inner<real_t>(va, va);
        }
        
        // Temperature at t(n).
        data.temperature = properties::temperature<T>(particles, data.ekin);
        
        // Time.
        data.t = counter * dt;
        
        return data;        
    }
    
    LeapFrog::LeapFrog(sim_param_ptr_t simulationParameters) :
        simulationParameters_{std::move(simulationParameters)} {
    }
    
    SimulationData 
    LeapFrog::displace(std::vector<std::shared_ptr<Particle>> &particles) const
    {
        static bool setup = false;
        static stime_t dt{0.0};
        static std::size_t counter = 0.0;
        
        counter += 1;
        
        if ( !setup ) {
            dt = param.get<real_t>("timestep");        
        }
        
        // Forces and energies.
        auto result = interactor_->interact(param, at);
        
        // Displace.
        SimulationData data = at->doWithAll<SimulationData>([] (const std::vector<atom_ptr_t>& atoms) {
            return displace_<Atom>(dt, atoms);
        });
        
        // Save simulation data.
        data.bepot = result.first;
        data.nbepot = result.second;
        data.t = counter * dt;
        
        return data;
    }
    
    std::string 
    LeapFrog<Atomistic>::id() const
    {
        return conf::LEAP_FROG;
    }
    
    LeapFrog<CoarseGrained>::LeapFrog(const cg_interactor_ptr_t& interactor) : 
        interactor_{interactor}
    {        
    }
    
    SimulationData 
    LeapFrog<CoarseGrained>::displace(const sim_param_t& param, 
                                      const cg_mod_ptr_t& cg) const
    {
        static bool setup = false;
        static stime_t dt{0.0};
        static std::size_t counter = 0.0;
        
        counter += 1;
        
        if ( !setup ) {
            dt = param.get<real_t>("timestep");        
        }
        
        // Forces and energies.
        auto result = interactor_->interact(param, cg);
        SimulationData data = cg->doWithAll<SimulationData>([] (const std::vector<bead_ptr_t>& beads) {
            return displace_<Bead>(dt, beads);
        });
        
        // Save simulation data.
        data.bepot = result.first;
        data.nbepot = result.second;
        data.t = counter * dt;
        
        return data;
    }
    
    std::string 
    LeapFrog<CoarseGrained>::id() const
    {
        return conf::LEAP_FROG;
    }
}

