/*
 * Author: André H. Juffer.
 * Created on 24/06/22, 13:34.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/s1-dpd.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <random>

namespace simploce {
    namespace s1_dpd {

        /**
         * Step 1: Displace velocities due to random and dissipative forces. See Nikunen et al., Table 4, p. 414.
         * @param pairLists Particle pair lists.
         * @param dt Timestep.
         * @param bc Boundary condition.
         * @param gamma Friction
         * @param cutoff Cutoff distance.
         * @param temperature Temperature.
         */
        static void
        randomDissipativeForces_(const PairLists &pairLists,
                                 const stime_t &dt,
                                 const bc_ptr_t &bc,
                                 real_t gamma,
                                 const dist_t &cutoff,
                                 const temperature_t &temperature) {
            util::Logger logger("simploce::s1_dpd::randomDissipativeForces_()");
            logger.trace("Entering.");

            std::random_device rd;
            std::mt19937_64 generator(rd());
            std::normal_distribution<real_t> normalDistribution{0.0, 1.0}; // Standard Wiener/Brownian

            static real_t KB = 1.0;
            static std::size_t counter{0};
            static real_t sigma = std::sqrt(2.0 * gamma * KB * temperature());
            static real_t sqrt_dt = std::sqrt(dt());
            static auto eps = real_t(std::numeric_limits<float>::min());

            if (counter == 0) {
                logger.info(std::to_string(sigma) + ": Strength (sigma) of random force.");
            }
            counter += 1;

            for (auto& pp: pairLists.particlePairList()) {
                auto &p1 = pp.first;
                auto &p2 = pp.second;
                auto &particle1 = *p1;
                auto &particle2 = *p2;

                auto r1 = particle1.position();    // Position at t(n+1).
                auto r2 = particle2.position();    // Position at t(n+1).

                auto rij = bc->apply(r1, r2);
                auto dist = norm<real_t>(rij);
                if (dist < cutoff()) {
                    // Random and dissipative pair forces.
                    auto mass1 = particle1.mass();
                    auto mass2 = particle2.mass();
                    auto v1 = particle1.velocity();    // Velocity at t(n+1).
                    auto v2 = particle2.velocity();
                    auto uv = (dist > eps ? (rij / dist) : util::randomUnit());
                    auto vij = v1 - v2;                  // Difference velocity.
                    auto ip = inner<real_t>(uv, vij);
                    auto w = 1.0 - dist / cutoff();               // (1-r/r_c)
                    auto w2 = w * w;                              // (1-r/r_c)^2
                    auto W = normalDistribution(generator);     // Random number.
                    auto t1 = -0.5 * gamma * w2 * ip * dt() + 0.5 * sigma * w * W * sqrt_dt;
                    auto t2 = -t1;
                    for (auto k = 0; k != 3; ++k) {
                        v1[k] += t1 / mass1() * uv[k];
                        v2[k] += t2 / mass2() * uv[k];
                    }
                    auto t3 = 0.5 * sigma * w * W * sqrt_dt;
                    t3 += -0.5 * (gamma * w2 * dt()) * (ip + sigma * w * W * sqrt_dt) / (1.0 + gamma * w2 * dt());
                    auto t4 = -t3;
                    for (auto k = 0; k != 3; ++k) {
                        v1[k] += t3 / mass1() * uv[k];
                        v2[k] += t4 / mass2() * uv[k];
                    }
                    particle1.velocity(v1);
                    particle2.velocity(v2);
                }
            }

            logger.trace("Leaving.");
        }

        /**
         * Steps 2 and 3: Displace both position and velocities. Forces on particle must be conservative forces only.
         * See Nikunen et al., Table 4, p. 414.
         * @param dt Time step.
         * @param particles Particles.
         */
        static void
        displace_(const stime_t &dt,
                  const std::vector<p_ptr_t> &particles) {
            util::Logger logger("simploce::s1_dpd::displace_()");
            logger.trace("Entering.");

            for (const auto& p : particles) {
                auto &particle = *p;
                auto mass = particle.mass();
                auto v = particle.velocity();
                auto r = particle.position();
                auto f = particle.force();
                for (auto k = 0; k != 3; ++k) {
                    v[k] += 0.5 * f[k] * dt() / mass();
                    r[k] += v[k] * dt();
                }
                particle.position(r);
                particle.velocity(v);
            }

            logger.trace("Leaving.");
        }

        /**
         * Final step 5
         * @param dt Time step.
         * @param particles Particles.
         * @return Simulation data.
         */
        static SimulationData
        displaceVelocity_(const stime_t &dt,
                          const std::vector<p_ptr_t> &particles) {
            util::Logger logger("simploce::s1_dpd::displaceVelocity_()");
            logger.trace("Entering.");

            SimulationData data;

            // Kinetic energy at t(n+1).
            data.kinetic = 0.0;
            for (const auto& p : particles) {
                auto &particle = *p;
                auto mass = particle.mass();
                auto v = particle.velocity();
                auto f = particle.force();
                for (auto k = 0; k != 3; ++k) {
                    v[k] += 0.5 * f[k] * dt() / mass();
                }
                particle.velocity(v);
                data.kinetic += 0.5 * mass() * inner<real_t>(v, v);
            }

            // Instantaneous temperature at t(n+1).
            data.temperature = properties::kineticTemperature(particles, data.kinetic, true);

            logger.trace("Leaving.");
            return data;
       }

    }

    S1_DPD::S1_DPD(param_ptr_t param,
                   interactor_ptr_t interactor,
                   bc_ptr_t bc,
                   units::dpd_ptr_t dpdUnits) :
            param_{std::move(param)}, interactor_{std::move(interactor)}, bc_{std::move(bc)},
            dpdUnits_{std::move(dpdUnits)} {
        util::Logger logger("simploce::S1_DPD::S1_DPD()");
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
    S1_DPD::displace(const p_system_ptr_t &particleSystem) const {
        static util::Logger logger("simploce::S1_DPD::displace()");
        logger.trace("Entering.");

        static std::size_t counter = 0;

        // Assuming DPD units
        static stime_t dt =  param_->get<real_t>("simulation.timestep");
        static temperature_t temperature = param_->get<real_t>("simulation.temperature");
        static auto gamma = param_->get<real_t>("simulation.displacer.dpd.gamma");
        static auto lambda = param_->get<real_t>("simulation.displacer.dpd.lambda");
        static dist_t cutoff = param_->get<real_t>("simulation.forces.cutoff");
        static auto isMesoscale = param_->get<bool>("simulation.mesoscale");
        static auto conservativeForces =
                param_->get<bool>("simulation.forces.conservative", true);
        if (!isMesoscale) {
            logger.warn("Simulation parameters may not be suitable for a DPD simulation.");
        }

        counter += 1;

        if ( counter == 1) {
            logger.info("Dissipative Particle Dynamics (S1).");
            logger.debug(util::toString(dt) + ": Time step.");
            logger.debug(util::toString(temperature) + ": Temperature.");
            logger.debug(util::toString(gamma) + ": gamma.");
            logger.debug(util::toString(lambda) + ": lambda.");
            logger.debug(util::toString(cutoff) + ": Cutoff distance.");

            // Initial conservative forces.
            interactor_->interact(particleSystem);
            if (!conservativeForces) {
                logger.warn("Conservative forces are ignored. Only thermostat is employed.");
                particleSystem->resetForces();
            }
        }

        // Steps 1 to 3: Displacement due to random and dissipative forces, followed by displacement
        // due to conservative forces.
        particleSystem->doWithDisplaceables<void>([this] (const std::vector<p_ptr_t>& particles) {
            auto pairLists = this->interactor_->pairLists();
            s1_dpd::randomDissipativeForces_(pairLists, dt, this->bc_, gamma, cutoff, temperature);
            s1_dpd::displace_(dt, particles);
        });

        // Step 4: Conservative forces and energies at t(n+1).
        auto result = this->interactor_->interact(particleSystem);
        if ( !conservativeForces ) {
            particleSystem->resetForces();
        }

        // Final step 5. Displace velocities.
        auto data =
                particleSystem->doWithDisplaceables<SimulationData>([] (
                        const std::vector<p_ptr_t>& particles) {
            auto data = s1_dpd::displaceVelocity_(dt, particles);
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

