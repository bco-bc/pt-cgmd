/*
 * Author: André H. Juffer.
 * Created on 17/11/2021, 15:55.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/simulation.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/displacer.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/simulation/s-util.hpp"

namespace simploce {

    Simulation::Simulation(param_ptr_t param,
                           p_system_ptr_t particleSystem,
                           displacer_ptr_t displacer,
                           bc_ptr_t bc) :
            param_{std::move(param)},
            particleSystem_{std::move(particleSystem)},
            displacer_{std::move(displacer)},
            bc_{std::move(bc)} {
    }

    void Simulation::perform(std::ofstream& trajectoryStream,
                             std::ofstream& dataStream) {
        static util::Logger logger("simploce::Simulation::perform()");
        logger.trace("Entering");

        // Set up.
        std::size_t numberAccepted = 0;
        auto box = particleSystem_->box();
        auto nSteps = param_->get<int>("simulation.nsteps");
        auto nWrite = param_->get<int>("simulation.nwrite");
        auto referenceTemperature = param_->get<real_t>("simulation.temperature");
        auto mesoscopic = param_->get<bool>("simulation.mesoscale");
        auto removeCenterOfMassMotion = param_->get<bool>("simulation.remove-com-motion", false);
        auto nRemoveCenterOfMassMotion =
                param_->get<int>("simulation.nremove-com-motion", nSteps + 1);
        auto scaleVelocities = param_->get<bool>("simulation.scale-velocities", false);
        auto nScaleVelocities =
                param_->get<int>("simulation.nscale-velocities", nSteps + 1);
        auto maxTemperatureDifference =
                param_->get<real_t>("simulation.relative-temperature-difference", 10.0);  // in %.

        logger.debug(std::to_string(nSteps) + ": Requested number of simulation steps.");
        logger.debug(std::to_string(nWrite) + ": Number of steps between writing simulation data.");
        logger.info(std::to_string(removeCenterOfMassMotion) + ": Remove center of mass motion?");
        logger.info(std::to_string(scaleVelocities) + ": Scale velocities?");

        if ( scaleVelocities ) {
            logger.info(std::to_string(referenceTemperature) +
                        ": Reference temperature for scaling velocities.");
            logger.debug(std::to_string(nScaleVelocities) + ": Number of steps between scaling velocities.");
            logger.debug(std::to_string(maxTemperatureDifference) +
                        ": Allowed relative difference (%) between actual and reference temperature.");
        }
        if ( removeCenterOfMassMotion ) {
            logger.debug(std::to_string(nRemoveCenterOfMassMotion) +
                        ": Number of steps between removing center of mass motion.");
        }

        // Perform simulation.
        logger.info("Simulation ongoing.");
        for (int counter = 1; counter <= nSteps; ++counter) {
            //logger.debug("Step #: " + std::to_string(counter));

            // Update positions and velocities.
            SimulationData data = displacer_->displace(particleSystem_);
            if ( data.accepted ) {
                numberAccepted += 1;
            }

            // Apply boundary conditions to velocities of displaceable particles (not frozen), if any.
            particleSystem_->doWithAll<void>([this] (const std::vector<p_ptr_t>& particles) {
                for (auto& p: particles) {
                    if (!p->frozen()) {
                        auto v = p->velocity();
                        auto r = p->position();
                        v = this->bc_->apply(v, r);
                        p->velocity(v);
                    }
                }
            });

            if ( counter % nWrite == 0) {
                data.pressure =
                        particleSystem_->doWithAll<pressure_t>([box, data] (
                                const std::vector<p_ptr_t>& particles) {
                    return properties::pressure(particles, data.temperature, box);
                });
                data.acceptanceRatio = real_t(numberAccepted)/real_t(counter) * 100.0;
                dataStream << std::setw(conf::INTEGER_WIDTH) << counter << conf::SPACE << data << std::endl;
                particleSystem_->writeState(trajectoryStream);
                trajectoryStream.flush();
                dataStream.flush();
            }
            if (removeCenterOfMassMotion) {
                if (counter % nRemoveCenterOfMassMotion == 0) {
                    logger.debug(std::to_string(counter) + ": Removing center of mass motion.");
                    util::removeCenterOfMassMotion(particleSystem_);
                }
            }
            if (scaleVelocities) {
                if (counter % nScaleVelocities == 0) {
                    auto difference =
                        std::fabs(referenceTemperature - data.temperature()) / referenceTemperature * 100.0;
                    if (difference > maxTemperatureDifference) {
                        logger.warn(std::to_string(counter) + ": Scaling velocities after displacement.");
                        util::scaleVelocities(particleSystem_, referenceTemperature, mesoscopic);
                    }
                }
            }
        }
        logger.info("Simulation completed.");

        logger.trace("Leaving.");
    }

}

