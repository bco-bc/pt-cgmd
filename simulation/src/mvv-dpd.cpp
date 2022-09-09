/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on April 28, 2022, 15:15 PM
 */

#include "simploce/simulation/mvv-dpd.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/units/units-dpd.hpp"
#include "simploce/util/util.hpp"
#include <utility>
#include <stdexcept>
#include <random>

namespace simploce {
    namespace mvv_dpd {

        static std::vector<force_t> fis_{};      // Forces at time t(n).
        static std::vector<velocity_t> vis_{};   // Velocities at time t(n).

        /**
         * Displaces position and velocities, the latter requires "correction" afterwards.
         * @param dt Time step.
         * @param lambda Value of lambda.
         * @param particles Particles.
         */
        static void
        displace_(const stime_t &dt,
                  real_t lambda,
                  const std::vector<p_ptr_t> &particles) {
            static util::Logger logger("simploce::mvv_dpd::displace_()");
            logger.trace("Entering");

            static auto dt2 = dt() * dt();
            for (const auto &p: particles) {
                auto &particle = *p;
                auto index = particle.index();

                auto ri = particle.position();  // Position at time t(n).
                auto vi = particle.velocity();  // Velocity at time t(n).
                auto mass = particle.mass();
                auto fi = particle.force();        // Force at time t(n).

                // Save velocity and force at t(n) for later use.
                vis_[index] = vi;
                fis_[index] = fi;

                auto a1 = dt2 / (2.0 * mass());
                auto a2 = lambda * dt() / mass();
                static position_t rf;                     // Position at time t(n+1),
                static velocity_t vf;                     // "Uncorrected" velocity at time t(n+1).
                for (std::size_t k = 0; k != 3; ++k) {
                    rf[k] = ri[k] + dt() * vi[k] + a1 * fi[k];
                    vf[k] = vi[k] + a2 * fi[k];
                }

                // Update particle t(n+1).
                particle.position(rf);
                particle.velocity(vf);
            }

            logger.trace("Leaving");
        }

        static SimulationData
        correctVelocity_(const stime_t &dt,
                         const std::vector<p_ptr_t> &particles) {
            util::Logger logger("simploce::mvv_dpd::correctVelocity_()");
            logger.trace("Entering");

            SimulationData data;

            // Kinetic energy at t(n+1).
            data.kinetic = 0.0;

            for (const auto &p: particles) {
                auto &particle = *p;
                auto index = particle.index();

                auto mass = particle.mass();
                auto vi = vis_[index];               // Velocity at time t(n).
                auto fi = fis_[index];               // Force at time t(n)
                auto ff = particle.force();             // Force at time t(n+1).

                auto a1 = dt() / (2.0 * mass());
                static velocity_t vf{};                         // Velocity at time t(n+1).
                for (auto k = 0; k != 3; ++k) {
                    vf[k] = vi[k] + a1 * (fi[k] + ff[k]);
                }
                particle.velocity(vf);

                // Kinetic energy at time t(n+1).
                data.kinetic += 0.5 * mass() * inner<real_t>(vf, vf);
            }

            // Instantaneous temperature at t(n+1).
            data.temperature = properties::kineticTemperature(particles, data.kinetic, true);

            // Done.
            logger.trace("Leaving");
            return data;
        }

        static void
        randomDissipativeForces_(const PairLists &pairLists,
                                 const stime_t &dt,
                                 const bc_ptr_t &bc,
                                 real_t gamma,
                                 real_t weightFactor,
                                 const dist_t &cutoff,
                                 const temperature_t &temperature) {
            static util::Logger logger("simploce::mvv_dpd::randomDissipativeForces_()");
            logger.trace("Entering.");

            std::random_device rd;
            std::mt19937_64 generator(rd());
            std::normal_distribution<real_t> normalDistribution{0.0, 1.0}; // Standard Wiener/Brownian

            // Unit of temperature such that kB*T = n, where n is an integer. See units-mvv_dpd.hpp.
            static real_t KB = 1.0;
            static std::size_t counter{0};
            static real_t sigma = std::sqrt(2.0 * gamma * KB * temperature());
            static real_t factor = 1.0 / std::sqrt(dt());
            static auto eps = real_t(std::numeric_limits<float>::min());

            if (counter == 0) {
                logger.info(std::to_string(sigma) + ": Strength (sigma) of random force.");
                logger.info(std::to_string(weightFactor) + ": Weight factor s in (1.0-r/r_c)^s.");
            }
            counter += 1;

            for (const auto &pp: pairLists.particlePairList()) {
                auto &p1 = pp.first;
                auto &p2 = pp.second;
                auto &particle1 = *p1;
                auto &particle2 = *p2;

                auto r1 = particle1.position();    // Position at t(n+1).
                auto v1 = particle1.velocity();    // "Uncorrected" velocity at t(n+1).
                auto r2 = particle2.position();    // Position at t(n+1).
                auto v2 = particle2.velocity();    // "Uncorrected" velocity at t(n+1).

                auto rij = bc->apply(r1, r2);
                auto dist = norm<real_t>(rij);
                if (dist < cutoff()) {
                    // Random and dissipative pair forces.
                    auto uv = (dist > eps ? (rij / dist) : util::randomUnit());
                    auto vij = v1 - v2;                  // Difference velocity.
                    auto ip = inner<real_t>(uv, vij);
                    auto w = 1.0 - dist / cutoff();                // (1-r/r_c)
                    auto wD = std::pow(w, weightFactor);      // (1-r/r_c)^s
                    auto wR = std::sqrt(wD);                     // wD = (wR)^2
                    auto W = normalDistribution(generator);      // Random number.
                    static force_t randomF_1;                              // Random force on particle 1
                    static force_t dissipativeF_1;                         // Dissipative force on particle 1.
                    for (auto k = 0; k != 3; ++k) {
                        randomF_1[k] = sigma * wR * W * factor * uv[k];  // Eq (8) in Groot and Warren, 1997.
                        dissipativeF_1[k] = -gamma * wD * ip * uv[k];    // First of Eq (4) in Groot and Warren, 1997.
                    }

                    // Update the forces acting on particles 1 and 2.
                    auto f1 = particle1.force();
                    auto f2 = particle2.force();
                    for (auto k = 0; k != 3; ++k) {
                        auto total_1 = randomF_1[k] + dissipativeF_1[k];
                        f1[k] += total_1;
                        f2[k] -= total_1;
                    }
                    particle1.force(f1);
                    particle2.force(f2);
                }
            }

            logger.trace("Leaving");
        }

    }

    MVV_DPD::MVV_DPD(param_ptr_t param, interactor_ptr_t interactor, bc_ptr_t bc, units::dpd_ptr_t dpdUnits) :
            param_{std::move(param)}, interactor_{std::move(interactor)}, bc_{std::move(bc)},
            dpdUnits_{std::move(dpdUnits)} {
        util::Logger logger("simploce::MVV_DPD::MVV_DPD()");
        if (!param_) {
            util::logAndThrow(logger, "Missing simulation parameters.");
        }
        if (!interactor_) {
            util::logAndThrow(logger, "Missing interactor.");
        }
        if ( !bc_) {
            util::logAndThrow(logger, "Missing boundary conditions.");
        }
        if ( !dpdUnits_ ) {
            util::logAndThrow(logger, "Missing dpd units.");
        }
    }

    SimulationData
    MVV_DPD::displace(const p_system_ptr_t& particleSystem) const
    {
        static util::Logger logger("simploce::MVV_DPD::displace()");
        logger.trace("Entering.");

        static std::size_t counter = 0;

        // Assuming DPD units
        static stime_t dt =  param_->get<real_t>("simulation.timestep");
        static temperature_t temperature = param_->get<real_t>("simulation.temperature");
        static auto gamma = param_->get<real_t>("simulation.displacer.dpd.gamma");
        static auto lambda = param_->get<real_t>("simulation.displacer.dpd.lambda");
        static dist_t cutoff = param_->get<real_t>("simulation.forces.cutoff");
        static auto weightFactor =
                param_->get<real_t>("simulation.displacer.dpd.weight-factor",
                                                        1.0);
        static auto isMesoscale = param_->get<bool>("simulation.mesoscale");
        static auto conservativeForces =
            param_->get<bool>("simulation.forces.conservative", true);
        if (!isMesoscale) {
            logger.warn("Simulation parameters may not be suitable for a MVV_DPD simulation.");
        }

        counter += 1;

        if (counter == 1) {
            logger.info("Dissipative Particle Dynamics (Modified Velocity Verlet).");
            logger.debug(util::toString(dt) + ": Time step.");
            logger.debug(util::toString(temperature) + ": Temperature.");
            logger.debug(util::toString(gamma) + ": gamma.");
            logger.debug(util::toString(lambda) + ": lambda.");
            logger.debug(util::toString(cutoff) + ": Cutoff distance.");
            logger.debug(std::to_string(weightFactor) + ": Weight factor (s) in (1-r/rc)^s.");

            particleSystem->doWithDisplaceables<void>([] (const std::vector<p_ptr_t>& particles) {
                auto numberOfParticles = particles.size();
                mvv_dpd::fis_ = std::vector<force_t>(numberOfParticles, force_t{});
                mvv_dpd::vis_ = std::vector<velocity_t>(numberOfParticles, velocity_t{});
            });

            // Initial conservative forces.
            interactor_->interact(particleSystem);

            if (!conservativeForces) {
                logger.warn("Conservative forces are ignored. Only thermostat is employed.");
                particleSystem->resetForces();
            }

            // Add random and dissipative forces.
            const auto& pairLists = interactor_->pairLists();
            mvv_dpd::randomDissipativeForces_(pairLists, dt, bc_, gamma, weightFactor, cutoff, temperature);
        }

        // Displace position and velocities, the latter are the "uncorrected" velocities.
        particleSystem->doWithDisplaceables<void>([](const std::vector<p_ptr_t>& particles) {
            mvv_dpd::displace_(dt, lambda, particles);
        });

        // Compute conservative forces and potential energy at t(n+1) using positions at t(n+1).
        auto result = interactor_->interact(particleSystem);
        if ( !conservativeForces ) {
            particleSystem->resetForces();
        }

        // Add random and dissipative forces using positions and "uncorrected" velocities at time t(n+1).
        const auto& pairLists = interactor_->pairLists();
        mvv_dpd::randomDissipativeForces_(pairLists, dt, bc_, gamma, weightFactor, cutoff, temperature);

        // "Correct" particle velocities.
        SimulationData data =
                particleSystem->doWithDisplaceables<SimulationData>([] (const std::vector<p_ptr_t>& particles) {
            auto data = mvv_dpd::correctVelocity_(dt, particles);
            data.totalMomentum = norm<real_t>(properties::linearMomentum(particles));
            return data;
        });

        // Save simulation data at t(n+1).
        data.bonded = std::get<0>(result);
        data.nonBonded = std::get<1>(result);
        data.external = std::get<2>(result);
        data.t = counter * dt;
        data.accepted = true;

        logger.trace("Leaving.");
        return data;
    }
}