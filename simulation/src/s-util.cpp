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

namespace simploce {
    namespace util {

        void scaleVelocities(const p_system_ptr_t& particleSystem,
                             const temperature_t& referenceTemperature,
                             bool isMesoscale) {
            static util::Logger logger("simploce::util::scaleVelocities()");
            logger.trace("Entering.");

            logger.debug(std::to_string(isMesoscale) + ": Mesoscale particle system?");
            logger.debug(std::to_string(referenceTemperature()) + ": Reference temperature.");

            particleSystem->doWithDisplaceables<void>([referenceTemperature, isMesoscale] (
                    const std::vector<p_ptr_t>& displaceables) {
                energy_t ekin{0.0};
                for (const auto& p: displaceables) {
                    auto v = p->velocity();
                    ekin += 0.5 * p->mass()() * inner<real_t>(v, v);
                }
                auto temperature =
                        properties::kineticTemperature(displaceables, ekin, isMesoscale);
                auto scaling = std::sqrt(referenceTemperature()/ temperature());
                logger.debug(std::to_string(temperature()) + ": Instantaneous (kinetic) temperature.");
                logger.debug(std::to_string(scaling) + ": Temperature scaling factor.");

                for (const auto& p: displaceables) {
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

