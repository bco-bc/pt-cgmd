/*
 * File:   leap-frog.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 5:12 PM
 */

#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include <vector>

namespace simploce {
    
    /*
     * @param T particle pointer type.
     * @param dt Time step. 
     * @param particles Particles.
     * @return Kinetic energy, temperature.
     */
    static SimulationData 
    displace_(const stime_t dt, 
              const std::vector<p_ptr_t>& particles)
    {
        static std::size_t counter = 0;
        
        // Assume current step n-1/2 at time t(n-1/2).
        
        // Compute linear momentum and position, plus kinetic energy.
        SimulationData data;
        for (auto ptr : particles) {
            auto &particle = *ptr;
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
            data.kinetic += 0.5 * mass() * inner<real_t>(va, va);
        }
        
        // Temperature at t(n).
        data.temperature = properties::temperature(particles, data.kinetic);
        
        return std::move(data);
    }
    
    LeapFrog::LeapFrog(sim_param_ptr_t simulationParameters,
                       interactor_ptr_t interactor) :
        simulationParameters_{std::move(simulationParameters)}, interactor_{interactor} {
    }
    
    SimulationData 
    LeapFrog::displace(const p_system_ptr_t& particleSystem) const
    {
        static stime_t dt = simulationParameters_->get<real_t>("simulation.timestep");
        static std::size_t counter = 0.0;
        
        counter += 1;
        
        // Forces and energies.
        auto result = interactor_->interact(particleSystem);
        
        // Displace.
        SimulationData data =
                particleSystem->doWithAll<SimulationData>([] (const std::vector<p_ptr_t>& all) {
            return std::move(displace_(dt, all));
        });
        
        // Save simulation data.
        data.bonded = result.first;
        data.nonBonded = result.second;
        data.t = counter * dt;
        
        return data;
    }
    
}

