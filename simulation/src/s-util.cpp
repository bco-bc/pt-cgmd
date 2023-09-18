/*
 * Author: André H. Juffer.
 * Created on 15/06/2022.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/s-util.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/units/units-mu.hpp"

namespace simploce {
    namespace util {

        dist_t
        computePairListCutoff(const param_ptr_t& param,
                              const p_system_ptr_t & particleSystem) {

            util::Logger logger("simploce::util::computePairListCutoff()");
            logger.trace("Entering.");

            temperature_t temperature = param->get<real_t>("simulation.temperature");
            auto mesoscopic = param->get<bool>("simulation.mesoscale", false);
            stime_t dt = param->get<real_t>("simulation.timestep", 0.001);
            auto cutoff = param->get<real_t>("simulation.forces.cutoffLR");
            auto nPairLists = param->get<int>("simulation.npairlists");

            // Compute typical displacement.
            auto numberOfParticles = particleSystem->numberOfParticles();
            auto KB = mesoscopic ? 1 : units::mu<real_t>::KB;
            mass_t mass = particleSystem->mass();
            mass /= real_t(numberOfParticles);
            energy_t eKin = KB * temperature() / mass();  // Equipartition, 1D.
            auto v = std::sqrt(2.0 * eKin() / mass());
            auto displacement = v * dt;

            // Verlet cutoff distance for pair lists.
            // See https://en.wikipedia.org/wiki/Verlet_list
            dist_t pairListCutoff = cutoff + 2.0 * nPairLists * displacement();

            logger.debug(std::to_string(temperature()) + ": Temperature.");
            logger.debug(std::to_string(dt()) + ": Time step.");
            logger.debug(std::to_string(eKin()) +
                          ": Estimated kinetic energy per particle at the given temperature.");
            logger.debug(std::to_string(mass()) + ": Average mass.");
            logger.debug(std::to_string(mesoscopic) + ": Mesoscopic simulation?");
            logger.debug(std::to_string(cutoff) + ": Cutoff distance for long-ranged forces.");
            logger.debug(std::to_string(v) + ": Average velocity.");
            logger.debug(std::to_string(displacement()) + ": Typical displacement in single step.");
            logger.debug(std::to_string(nPairLists) + ": Number of steps between update pair lists.");
            logger.info(std::to_string(pairListCutoff()) + ": Cutoff distance for pair lists generation.");

            logger.trace("Leaving.");
            return pairListCutoff;
        }

        void scaleVelocities(const p_system_ptr_t& particleSystem,
                             const temperature_t& referenceTemperature,
                             bool mesoscopic) {
            static util::Logger logger("simploce::util::scaleVelocities()");
            logger.trace("Entering.");

            logger.debug(std::to_string(mesoscopic) + ": Mesoscopic particle system?");
            logger.debug(std::to_string(referenceTemperature()) + ": Reference temperature.");

            particleSystem->doWithAll<void>([referenceTemperature, mesoscopic] (
                    const std::vector<p_ptr_t>& all) {
                energy_t ekin{0.0};
                for (const auto& p: all) {
                    auto v = p->velocity();
                    ekin += 0.5 * p->mass()() * inner<real_t>(v, v);
                }
                auto temperature =
                        properties::kineticTemperature(all, ekin, mesoscopic);
                auto scaling = std::sqrt(referenceTemperature()/ temperature());
                logger.debug(std::to_string(temperature()) + ": Instantaneous (kinetic) temperature.");
                logger.debug(std::to_string(scaling) + ": Temperature scaling factor.");

                for (const auto& p: all) {
                    auto v = p->velocity();
                    v *= scaling;
                    p->velocity(v);
                }
            });

            logger.trace("Leaving.");
        }

        void removeCenterOfMassMotion(const p_system_ptr_t& particleSystem) {
            util::removeOverallLinearMomentum(particleSystem);
        }
    }
}

