/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 5:39 PM
 */

#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <utility>
#include <vector>

namespace simploce {
    namespace vv {

        /*
         * Displaces particle positions.
         * @param dt Time step.
         * @param particles Particles.
         * @return Forces at time t(n).
         */
        static std::vector<force_t>
        displace_(const stime_t &dt,
                  const std::vector<p_ptr_t> &particles) {
            // Initial forces, forces at time t(n).
            std::vector<force_t> fis(particles.size(), force_t{0.0, 0.0, 0.0});

            // Displace particles: Positions.
            for (std::size_t index = 0; index != particles.size(); ++index) {
                Particle &particle = *particles[index];

                mass_t mass = particle.mass();
                force_t fi = particle.force();                 // Force (kJ/(mol nm) = (u nm)/(ps^2))
                                                               // at time t(n).
                position_t ri = particle.position();           // Position at time t(n).
                velocity_t vi = particle.velocity();           // Velocity at time t(n).

                // Save force at t(n) for later use.
                fis[index] = fi;

                real_t a1 = dt() / (2.0 * mass());
                real_t a2 = dt() * a1;
                static position_t rf{};
                for (std::size_t k = 0; k != 3; ++k) {
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
        displaceVelocity_(const stime_t &dt,
                          const std::vector<force_t> &fis,
                          const std::vector<p_ptr_t> &particles) {
            SimulationData data;

            // Kinetic energy at t(n+1).
            data.kinetic = 0.0;

            // Displace particles: Momenta/velocities.
            for (std::size_t index = 0; index != particles.size(); ++index) {
                Particle &particle = *particles[index];

                mass_t mass = particle.mass();
                force_t fi = fis[index];                       // Force (kJ/(mol nm) = (u nm)/(ps^2))
                // at time t(n).
                force_t ff = particle.force();                 // Force (kJ/(mol nm) = (u nm)/(ps^2))
                // at time t(n+1).

                real_t a1 = dt() / (2.0 * mass());
                velocity_t vi = particle.velocity();           // Velocity (nm/ps) at time t(n).
                static velocity_t vf{};                        // Velocity at time t(n+1).
                for (std::size_t k = 0; k != 3; ++k) {
                    vf[k] = vi[k] + a1 * (fi[k] + ff[k]);
                }

                // Save velocity
                particle.velocity(vf);

                // Kinetic energy at t(n+1).
                data.kinetic += 0.5 * mass() * inner<real_t>(vf, vf);
            }

            // Instantaneous temperature at t(n+1).
            data.temperature = properties::kineticTemperature(particles, data.kinetic);

            return data;
        }
    }
       
    VelocityVerlet::VelocityVerlet(param_ptr_t param,
                                   interactor_ptr_t interactor) :
            param_{std::move(param)}, interactor_{std::move(interactor)} {
        if (!param_) {
            throw std::domain_error("VelocityVerlet: Missing simulation parameters.");
        }
        if (!interactor_) {
            throw std::domain_error("VelocityVerlet: Missing interactor.");
        }
    }
        
    SimulationData 
    VelocityVerlet::displace(const p_system_ptr_t& particleSystem) const
    {
        static util::Logger logger("simploce::VelocityVerlet::displace(const p_system_ptr_t& particleSystem)");
        logger.trace("Entering.");

        static std::size_t counter = 0;
        static stime_t dt =  param_->get<real_t>("simulation.timestep");
        static std::vector<force_t> fis{};

        if (counter == 0) {
            logger.debug("Time step: " + std::to_string(dt()));
        }
        
        counter += 1;

        if ( counter == 1) {
            interactor_->interact(particleSystem);  // Initial forces at t(n).
        }
        
        // Displace particle positions.
        particleSystem->doWithAll<void>([] (const std::vector<p_ptr_t>& particles) {
            fis = vv::displace_(dt, particles);
        });
        
        // Compute forces and potential energy at t(n+1) using positions at t(n+1).
        auto result = interactor_->interact(particleSystem);
        
        // Displace particle momenta.
        SimulationData data = particleSystem->doWithAll<SimulationData>([] (const std::vector<p_ptr_t>& particles) {
            auto data = vv::displaceVelocity_(dt, fis, particles);
            data.totalMomentum = norm<real_t>(properties::linearMomentum(particles));
            return data;
        });
        
        // Save simulation data at t(n+1).
        data.bonded = std::get<0>(result);
        data.nonBonded = std::get<1>(result);
        data.external = std::get<2>(result);
        data.t = counter * dt;
        data.accepted = true;

        logger.trace("Leaving");
        return data;
    }
}
