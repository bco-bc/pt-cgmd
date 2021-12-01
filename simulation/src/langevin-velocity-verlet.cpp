/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 9:47 AM
 */

#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <random>
#include <cmath>
#include <array>
#include <cassert>

namespace simploce {
    
    // Helpers.
    static std::vector<real_t> FC_{};
    static std::vector<real_t> B_{};
    static std::vector<real_t> A1_{};
    static std::vector<real_t> A2_{};
    static std::vector<real_t> strengths_{};
    
    static std::vector<force_t> fis_{};      // Forces at time t(n).
    static std::vector<position_t> ris_{};   // Positions at time t(n).
    
    // Random vector W, each element is an array of size 3.
    static std::random_device rd_{};
    static std::mt19937 gen_(rd_());
    static std::normal_distribution<real_t> dis_{0.0, 1.0}; // Standard Wiener/Brownian
    static std::vector<std::array<real_t, 3>> W_{};
    
    /**
     * Calculation of helper values.
     * @param dt Time step.
     * @param particles Particles.
     */
    static void
    setupHelpers_(const stime_t& dt,
                  const temperature_t& temperature,
                  real_t gamma,
                  const std::vector<p_ptr_t>& particles)
    {
        static util::Logger logger("simploce::LangevinVelocityVerlet::setupHelpers_");

        static real_t kT = units::mu<real_t>::KB * temperature();
        std::size_t nParticles = particles.size();
        
        FC_ = std::vector<real_t>(nParticles, 0.0);
        B_ = std::vector<real_t>(nParticles, 0.0);
        A1_ = std::vector<real_t>(nParticles, 0.0);
        A2_ = std::vector<real_t>(nParticles, 0.0);
        strengths_ = std::vector<real_t>(nParticles, 0.0);
        
        fis_ = std::vector<force_t>(nParticles, force_t{});
        ris_ = std::vector<position_t>(nParticles, position_t{});
        
        W_ = std::vector<std::array<real_t, 3>>(nParticles, std::array<real_t, 3>{0.0, 0.0, 0.0});
        auto value = util::seedValue<std::size_t>();
        gen_.seed(value);
    
        for (const auto& p : particles) {
            auto& particle = *p;
            
            auto index = particle.index();
            
            
            mass_t mass = particle.mass();                        // In u.
            real_t fc =  mass() * gamma;                          // Friction 
            FC_[index] = fc;                                      // coefficient in u/ps.
                                                                  
            // All kinds of constant factors for each particle.
            real_t a1 = dt() / (2.0 * mass());                    // ps/u
            A1_[index] = a1;
            real_t a2 = a1 * dt();                                // ps^2/u
            A2_[index] = a2;
            real_t strength = std::sqrt( dt() * 2.0 * fc * kT );  // In (u nm) / ps.
            strengths_[index] = strength;
            real_t b = 1.0 / ( 1.0 + fc * a1);                    // No units.
            B_[index] = b;

            // Validate.
            real_t f1 = fc * a1;
            if ( f1 > 1) {
                logger.warn("Factor = fc * dt / (2.0 * m) = " + util::toString(f1) +
                            " > 1 may represent an unphysical regime.");
            }
        }                
    }

    /**
     * Displace particle position.
     * @param dt Time step.
     * @param particles Particles.
     */
    static void
    displacePosition_(const stime_t& dt,
                      const std::vector<p_ptr_t>& particles)
    {        
        for (const auto& p : particles) {
            auto& particle = *p;
            
            auto index = particle.index();
    
            // Update position, not velocity.
            static std::array<real_t, 3> w;            // Random vector at t(n+1).
            w[0] = dis_(gen_);
            w[1] = dis_(gen_);
            w[2] = dis_(gen_);
            W_[index] = w;                             // Save for updating velocities.
        
            force_t fi = particle.force();             // Force (kJ/(mol nm) = 
            fis_[index] = fi;                          // (u nm)/(ps^2)) at time t(n)
      
            velocity_t vi = particle.velocity();       // Velocity (nm/ps) at time t(n).
            position_t ri = particle.position();       // Position at time t(n).
            ris_[index] = ri;                          // Save for velocity update.
            static position_t rf;
            real_t b = B_[index];                      // No units.
            real_t a1 = A1_[index];                    // ps/u
            real_t a2 = A2_[index];                    // ps^2/u
            real_t strength = strengths_[index];
            for (std::size_t k = 0; k != 3; ++k) { 
                rf[k] = ri[k] +
                    b * dt() * vi[k] +
                    b * a2 * fi[k] +
                    b * a1 * strength * w[k];          // Position at time t(n+1).
            }
        
            // Update position.
            particle.position(rf);        
        }
    }

    static SimulationData
    displaceVelocity_(const std::vector<p_ptr_t>& particles) {
        SimulationData data;
        
        // Kinetic energy at t(n+1).
        data.kinetic = 0.0;
        
        // Displace particles: momenta.
        for (auto& p : particles) {
            auto& particle = *p;
            
            auto index = particle.index();
            
            mass_t mass = particle.mass();             // In u.
      
            std::array<real_t, 3> w = W_[index];       // Random vector at t(n+1).
            force_t fi = fis_[index];                  // Force (kJ/(mol nm) = (u nm)/(ps^2)) 
                                                       // at time t(n).
            position_t ri = ris_[index];               // Position at time t(n).
      
            force_t ff = particle.force();             // Force (kJ/(mol nm) = (u nm)/(ps^2))
                                                       // at time t(n+1).
            velocity_t vi = particle.velocity();       // velocity (nm/ps) at time t(n).
            position_t rf = particle.position();       // Position at time t(n+1).
      
            real_t fc = FC_[index];
            real_t a1 = A1_[index];
            real_t strength = strengths_[index];
      
            static velocity_t vf{};
            for (std::size_t k = 0; k != 3; ++k) {
                vf[k] = vi[k] +
                    a1 * ( fi[k] + ff[k] ) -
                    fc * ( rf[k] - ri[k] ) / mass() +
                    strength * w[k] / mass();          // Velocity at time t(n+1).
            }
      
            // Update velocity.
            particle.velocity(vf);
      
            // Kinetic energy at time t(n+1)..
            data.kinetic += 0.5 * mass() * inner<real_t>(vf, vf);
        }
        
        // Instantaneous temperature at t(n+1).
        data.temperature = properties::temperature(particles, data.kinetic);

        // Done.
        return data;        
        
    }
    
    LangevinVelocityVerlet::LangevinVelocityVerlet(sim_param_ptr_t simulationParameters,
                                                   interactor_ptr_t interactor) :
        simulationParameters_(std::move(simulationParameters)), interactor_{std::move(interactor)} {
    }
        
    SimulationData 
    LangevinVelocityVerlet::displace(const p_system_ptr_t& particleSystem) const
    {
        static stime_t dt = simulationParameters_->get<real_t>("simulation.timestep");
        static temperature_t temperature = simulationParameters_->get<real_t>("simulation.temperature", 298.15);
        static real_t gamma = simulationParameters_->get<real_t>("simulation.gamma", 0.5);
        static std::size_t counter = 0;

        counter += 1;

        if ( counter == 1 ) {
            particleSystem->doWithAll<void>([] (const std::vector<p_ptr_t>& atoms) {
                setupHelpers_(dt, temperature, gamma, atoms);
            });
            interactor_->interact(particleSystem);  // Initial forces.
        }
                
        // Displace atom positions.
        particleSystem->doWithAll<void>([] (const std::vector<p_ptr_t>& all) {
            displacePosition_(dt, all);
        });
        
        // Compute forces and potential energy at t(n+1) using positions at t(n+1).
        auto result = interactor_->interact(particleSystem);
        
        // Displace atom velocities.
        SimulationData data = particleSystem->doWithAll<SimulationData>([] (const std::vector<p_ptr_t>& all) {
            auto data = displaceVelocity_(all);
            data.totalMomentum = norm<real_t>(properties::linearMomentum(all));
            return std::move(data);
        });
        
        // Save simulation data
        data.bonded = std::get<0>(result);
        data.nonBonded = std::get<1>(result);
        data.external = std::get<2>(result);
        data.t = counter * dt;
        data.accepted = true;
        
        return data;
    }    

}

