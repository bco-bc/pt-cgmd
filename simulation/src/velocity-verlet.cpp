/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 5:39 PM
 */

#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include <utility>
#include <vector>

namespace simploce {
    
    /*
     * Displaces particle positions.
     * @param dt Time step.
     * @param particles Particles.
     * @return Forces at time t(n).
     */
    static std::vector<force_t>
    displacePosition_(const stime_t& dt,
                      const std::vector<p_ptr_t>& particles)
    {
        // Initial forces (forces at time t(n).
        std::vector<force_t> fis(particles.size(), force_t{});
        
        // Displace particles: Positions.
        for (std::size_t index = 0; index != particles.size(); ++index) {
            Particle& particle = *particles[index];

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
     * @param dt Time step.
     * @param fis Forces at t(n)
     * @param particles Particles.
     * @return Kinetic, potential energy, and temperature.
     */
    static SimulationData
    displaceMomentum_(const stime_t& dt,
                      const std::vector<force_t>& fis,
                      const std::vector<p_ptr_t>& particles)
    {
        SimulationData data;
        
        // Kinetic energy at t(n+1).
        data.kinetic = 0.0;
        
        // Displace particles: Momenta/velocities.
        for (std::size_t index = 0; index != particles.size(); ++index) {
            Particle& particle = *particles[index];
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
            data.kinetic += 0.5 * mass() * inner<real_t>(vf, vf);
        }
    
        // Instantaneous temperature at t(n+1).
        data.temperature = properties::temperature(particles, data.kinetic);
        
        return std::move(data);
    }
       
    VelocityVerlet::VelocityVerlet(sim_param_ptr_t simulationParameters,
                                   interactor_ptr_t interactor) :
        simulationParameters_{std::move(simulationParameters)}, interactor_{std::move(interactor)} {
    }
        
    SimulationData 
    VelocityVerlet::displace(const p_system_ptr_t& particleSystem) const
    {        
        static std::size_t counter = 0;
        
        static stime_t dt =  simulationParameters_->get<real_t>("timestep");
        static std::vector<force_t> fis{};
        
        counter += 1;
        
        // Displace atom positions.
        particleSystem->doWithAll<void>([] (const std::vector<p_ptr_t>& atoms) {
            fis = displacePosition_(dt, atoms);
        });
        
        // Compute forces and potential energy at t(n+1) using positions at t(n+1).
        auto result = interactor_->interact(particleSystem);
        
        // Displace atom momenta.
        SimulationData data = particleSystem->doWithAll<SimulationData>([] (const std::vector<p_ptr_t>& atoms) {
            return std::move(displaceMomentum_(dt, fis, atoms));
        });
        
        // Save simulation data.
        data.bonded = result.first;
        data.nonBonded = result.second;
        data.t = counter * dt;
        
        return data;
    }
}
