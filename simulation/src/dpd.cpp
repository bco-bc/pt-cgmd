/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on April 28, 2022, 15:15 PM
 */

#include "simploce/simulation/dpd.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/util.hpp"
#include <utility>
#include <stdexcept>
#include <random>

namespace simploce {

    static std::vector<force_t> fis_{};      // Forces at time t(n).
    static std::vector<velocity_t> vis_{};   // Velocities at time t(n).
    static std::random_device rd_{};
    static std::mt19937 gen_(rd_());
    static std::normal_distribution<real_t> dis_{0.0, 1.0}; // Standard Wiener/Brownian

    static void
    setupHelpers_(const stime_t& dt,
                  const temperature_t& temperature,
                  const std::vector<p_ptr_t>& particles)
    {
        std::size_t numberOfParticles = particles.size();
        fis_ = std::vector<force_t>(numberOfParticles, force_t{});
        vis_ = std::vector<velocity_t>(numberOfParticles, velocity_t{});
        auto value = util::seedValue<std::size_t>();
        gen_.seed(value);
    }

    /**
     * Displaces position.
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

           position_t ri = particle.position();  // Position at time t(n).
           velocity_t vi = particle.velocity();  // Velocity at time t(n).
           force_t fi = particle.force();        // Force at time t(n)
           fis_[index] = fi;

           static position_t rf;                 // Position at time t(n+1),
           for (std::size_t k = 0; k !=3; ++k) {
               rf[k] = ri[k] + dt() * vi[k] + 0.5 * dt() * dt() + fi[k];
           }

           // Update position.
           particle.position(rf);
       }
    }

    static void
    displaceVelocityUncorrected_(const stime_t& dt,
                                 real_t lambda,
                                 const std::vector<p_ptr_t>& particles) {
        for (const auto& p: particles) {
            auto& particle = *p;

            auto index = particle.index();

            mass_t mass = particle.mass();
            velocity_t vi = particle.velocity();   // Velocity at time t(n)
            vis_[index] = vi;
            force_t fi = fis_[index];              // Force at time t(n)

            static velocity_t vf;                  // Uncorrected velocity at time t(n+1).
            for (std::size_t k = 0; k !=3; ++k) {
                vf[k] = vi[k] + lambda * dt() * fi[k] / mass();   // Uncorrected velocity at time t(n+1)
            }

            // Store uncorrected velocity.
            particle.velocity(vf);
        }
    }

    static SimulationData
    correctDisplacedVelocity_(const stime_t& dt,
                              const std::vector<p_ptr_t>& particles)
    {
        SimulationData data;

        // Kinetic energy at t(n+1).
        data.kinetic = 0.0;

        for (const auto& p: particles) {
            auto& particle = *p;

            auto index = particle.index();
            mass_t mass = particle.mass();
            velocity_t vi = vis_[index];               // Velocity (nm/ps) at time t(n).
            force_t fi = fis_[index];                  // Force (kJ/(mol nm) = (u nm)/(ps^2))
                                                       // at time t(n)
            force_t ff = particle.force();             // Force (kJ/(mol nm) = (u nm)/(ps^2))
                                                       // at time t(n+1).
            static velocity_t vf{};
            for (std::size_t k = 0; k !=3; ++k) {
                 vf[k] = vi[k] + 0.5 * dt() * (fi[k] + ff[k]) / mass();
            }
            particle.velocity(vf);

            // Kinetic energy at time t(n+1).
            data.kinetic += 0.5 * mass() * inner<real_t>(vf, vf);

            // Done.
            return data;
        }
    }

    static void
    randomDissipativeForces_(const PairLists& pairLists,
                             const stime_t& dt,
                             const bc_ptr_t& bc,
                             real_t gamma,
                             const dist_t& cutoff,
                             const temperature_t& temperature)
    {
        static real_t sigma = std::sqrt(2.0 * gamma * units::mu<real_t>::KB * temperature());
        static real_t factor = 1.0 / std::sqrt(dt());

        for ( const auto& pp: pairLists.particlePairList() ) {
            auto& p1 = pp.first;
            auto& p2 = pp.second;

            auto& particle1 = *p1;
            auto& particle2 = *p2;

            auto r1 = particle1.position();
            auto v1 = particle1.velocity();
            auto r2 = particle2.position();
            auto v2 = particle2.velocity();
            auto rij = bc->apply(r1, r2);
            auto dist = norm<real_t>(rij);
            auto uv = rij / dist;    // Unit vector.
            auto v = v1 - v2;

            // Random and dissipative pair forces.
            auto ip = inner<real_t>(uv, v);
            auto  t = 1.0 - dist/cutoff();  // Should be >= 0.0
            t = t < 0.0 ? 0.0 : t;
            auto wD = t * t;
            auto wR = std::sqrt(wD);
            auto W = dis_(gen_);         // Random number drawn from normal distribution.
            static force_t randomF1;               // Random force on particle 1
            static force_t dissipativeF1;          // Dissipative force on particle 1.
            for (std::size_t k = 0; k != 3; ++k) {
                randomF1[k] = sigma * wR * W * factor * uv[k];
                dissipativeF1[k] = -gamma * wD * ip * uv[k];
            }

            // Update the forces acting on particle 1 and 2.
            auto f1 = particle1.force();
            auto f2 = particle2.force();
            for (std::size_t k = 1; k != 3; ++k) {
                auto total = randomF1[k] + dissipativeF1[k];
                f1[k] += total;
                f2[k] -= total;
            }
            particle1.force(f1);
            particle2.force(f2);
        }
    }

    DPD::DPD(param_ptr_t param, interactor_ptr_t interactor, bc_ptr_t bc) :
        param_{std::move(param)}, interactor_{std::move(interactor)}, bc_{std::move(bc)} {
        if (!param_) {
            throw std::domain_error("DPD: Missing simulation parameters.");
        }
        if (!interactor_) {
            throw std::domain_error("DPD: Missing interactor.");
        }
        if ( !bc_) {
            throw std::domain_error("DPD: Missing boundary conditions.");
        }
    }

    SimulationData
    DPD::displace(const p_system_ptr_t& particleSystem) const
    {
        static std::size_t counter = 0;
        static stime_t dt =  param_->get<real_t>("simulation.timestep");
        static temperature_t temperature = param_->get<real_t>("simulation.temperature");
        static auto gamma = param_->get<real_t>("simulation.gamma");
        static auto lambda = param_->get<real_t>("simulation.dpd.lambda");
        static dist_t cutoff = param_->get<real_t>("forces.nb.cutoff");

        counter += 1;

        if (counter == 1) {
            particleSystem->doWithDisplaceables<void>([] (const std::vector<p_ptr_t>& particles) {
                setupHelpers_(dt, temperature, particles);
            });
            interactor_->interact(particleSystem);  // Initial forces.

            // Add random and dissipative forces.
            const auto& pairLists = interactor_->pairLists();
            randomDissipativeForces_(pairLists, dt, bc_, gamma, cutoff, temperature);
        }

        // Displace position and velocities, the latter are the "uncorrected" velocities.
        particleSystem->doWithDisplaceables<void>([](const std::vector<p_ptr_t>& particles) {
            displacePosition_(dt, particles);
            displaceVelocityUncorrected_(dt, lambda, particles);
        });

        // Compute forces and potential energy at t(n+1) using positions at t(n+1).
        auto result = interactor_->interact(particleSystem);

        // Add random and dissipative forces using positions and uncorrected velocities at time t(n+1).
        const auto& pairLists = interactor_->pairLists();
        randomDissipativeForces_(pairLists, dt, bc_, gamma, cutoff, temperature);

        // Displace particle velocities.
        SimulationData data = particleSystem->doWithDisplaceables<SimulationData>([] (const std::vector<p_ptr_t>& particles) {
            auto data = correctDisplacedVelocity_(dt, particles);
            data.totalMomentum = norm<real_t>(properties::linearMomentum(particles));
            return data;
        });

        return data;
    }
}