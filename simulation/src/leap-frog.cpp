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
#include <utility>
#include <vector>
#include <stdexcept>

namespace simploce {
    
    /*
     * @param dt Time step. 
     * @param particles Particles.
     * @return Kinetic energy, temperature.
     */
    static SimulationData 
    displace_(const stime_t dt, 
              const std::vector<p_ptr_t>& particles)
    {
        // Assume current step n-1/2 at time t(n-1/2).
        
        // Compute linear momentum and position, plus kinetic energy.
        SimulationData data;
        for (const auto& ptr : particles) {
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
        data.temperature = properties::kineticTemperature(particles, data.kinetic);
        
        return data;
    }
    
    LeapFrog::LeapFrog(param_ptr_t param,
                       interactor_ptr_t interactor) :
            param_{std::move(param)}, interactor_{std::move(interactor)} {
        if (!param_) {
            throw std::domain_error("LeapFrog: Missing simulation parameters.");
        }
        if (!interactor_) {
            throw std::domain_error("LeapFrog: Missing interactor.");
        }
    }
    
    SimulationData 
    LeapFrog::displace(const p_system_ptr_t& particleSystem) const
    {
        static stime_t dt = param_->get<real_t>("simulation.timestep");
        static std::size_t counter = 0;
        
        counter += 1;
        
        // Forces and energies.
        auto result = interactor_->interact(particleSystem);
        
        // Displace.
        SimulationData data =
                particleSystem->doWithAll<SimulationData>([] (const std::vector<p_ptr_t>& particles) {
            auto data = displace_(dt, particles);
            data.totalMomentum = norm<real_t>(properties::linearMomentum(particles));
            return std::move(data);
        });

        
        // Save simulation data.
        data.bonded = std::get<0>(result);
        data.nonBonded = std::get<1>(result);
        data.external = std::get<2>(result);
        data.t = counter * dt;
        data.accepted = true;
        
        return data;
    }
    
}

