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
#include "simploce/simulation/interactor.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/simulation/s-util.hpp"

namespace simploce {

    Simulation::Simulation(param_ptr_t param,
                           p_system_ptr_t particleSystem,
                           spec_catalog_ptr_t catalog,
                           displacer_ptr_t displacer,
                           bc_ptr_t bc,
                           interactor_ptr_t interactor) :
            param_{std::move(param)},
            particleSystem_{std::move(particleSystem)},
            catalog_{std::move(catalog)},
            displacer_{std::move(displacer)},
            bc_{std::move(bc)},
            interactor_{std::move(interactor)} {
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
        auto removeCenterOfMassMotionAfter =
                param_->get<int>("simulation.nremove-com-motion", nWrite);
        auto scaleVelocities = param_->get<bool>("simulation.scale-velocities", false);
        auto nScaleVelocities =
                param_->get<int>("simulation.nscale-velocities", nSteps + 1);
        auto applyBCToVelocities = param_->get<bool>("simulation.bc.velocities.include", false);
        auto applyBcForVelocitiesToParticleGroups =
                param_->get<bool>("simulation.bc.velocities.particle-groups", false);
        auto maxTemperatureDifference =
                param_->get<real_t>("simulation.relative-temperature-difference", 10.0);  // in %.
        auto freezeBoundaryParticles = param_->get<bool>("simulation.freeze-boundary");
        auto freeze = param_->get<std::string>("simulation.freeze", std::string{});
        auto applyBCToPositions =  param_->get<bool>("simulation.bc.positions.include", false);

        logger.info(std::to_string(nSteps) + ": Requested number of simulation steps.");
        logger.info(std::to_string(nWrite) + ": Number of steps between writing simulation data.");
        logger.info(std::to_string(removeCenterOfMassMotion) + ": Remove center of mass motion?");
        logger.info(std::to_string(scaleVelocities) + ": Scale velocities?");
        logger.info(std::to_string(freezeBoundaryParticles) + ": Freeze boundary particles?");
        if (!freeze.empty())
            logger.info(freeze + ": Freeze particles of this specification.");
        logger.info(std::to_string(applyBCToVelocities) + ": Apply boundary condition to velocities?");
        logger.info(std::to_string(applyBCToPositions) + ": Apply boundary condition to positions?");

        if ( scaleVelocities ) {
            logger.info(std::to_string(referenceTemperature) +
                        ": Reference temperature for scaling velocities.");
            logger.debug(std::to_string(nScaleVelocities) + ": Number of steps between scaling velocities.");
            logger.debug(std::to_string(maxTemperatureDifference) +
                        ": Allowed relative difference (%) between actual and reference temperature.");
        }
        if ( removeCenterOfMassMotion ) {
            logger.debug(std::to_string(removeCenterOfMassMotionAfter) +
                         ": Number of steps passed before removing center of mass motion.");
        }
        if (freezeBoundaryParticles) {
            auto spec = catalog_->staticBP();
            particleSystem_->freeze(spec);
            logger.debug(spec->name() + ": Particles of this specification are frozen.");
            logger.info(std::to_string(particleSystem_->numberOfFrozenParticles()) + ": Number of frozen particles.");
        }
        if (!freeze.empty()) {
            auto spec = catalog_->lookup(freeze);
            particleSystem_->freeze(spec);
            logger.info(spec->name() + ": Particles of this specification are frozen.");
            logger.info(std::to_string(particleSystem_->numberOfFrozenParticles()) + ": Number of frozen particles.");
        }
        if (applyBCToVelocities) {
            if (applyBcForVelocitiesToParticleGroups) {
                logger.info(std::to_string(applyBcForVelocitiesToParticleGroups) +
                            ": Apply boundary conditions for velocities to particle groups.");
            }
        }

        // Perform simulation.
        logger.info("Simulation ongoing.");
        interactor_->initiate(particleSystem_);
        for (int counter = 1; counter <= nSteps; ++counter) {

            // Update positions and velocities.
            SimulationData data = displacer_->displace(particleSystem_);
            if ( data.accepted ) {
                numberAccepted += 1;
            }

            // Apply boundary conditions to velocities and/or positions of displaceable particles
            // (not frozen), if any.
            if (applyBCToVelocities) {
                // Apply to velocities.
                particleSystem_->doWithAllFreeGroups<void>([this, counter, applyBcForVelocitiesToParticleGroups](
                    const std::vector<p_ptr_t> &all,
                    const std::vector<p_ptr_t> &free,
                    const std::vector<pg_ptr_t> &groups) {
                    if (applyBcForVelocitiesToParticleGroups) {
                        if (counter == 1) {
                            logger.debug("Applying to BC for velocities to all particle groups.");
                        }
                        for (auto &g: groups) {
                            this->bc_->applyToVelocities(g);
                        }
                        for (auto &p: free) {
                            // Exclude frozen particles.
                            if (!p->frozen()) {
                                auto v = p->velocity();
                                auto r = p->position();
                                v = this->bc_->apply(v, r);
                                p->velocity(v);
                            }
                        }
                    } else {
                        if (counter == 1) {
                            logger.debug("Applying boundary condition to particle velocities.");
                        }
                        for (auto &p: all) {
                            // Exclude frozen particles.
                            if (!p->frozen()) {
                                auto v = p->velocity();
                                auto r = p->position();
                                v = this->bc_->apply(v, r);
                                p->velocity(v);
                            }
                        }
                    }
                });
            }

            if (applyBCToPositions) {
                // Apply to positions.
                particleSystem_->doWithAllFreeGroups<void>([this, counter](
                    const std::vector<p_ptr_t> &all,
                    const std::vector<p_ptr_t> &free,
                    const std::vector<pg_ptr_t> &groups) {
                    if (counter == 1) {
                        logger.debug("Applying boundary condition to particle positions.");
                    }
                    for (auto& p: all) {
                        if (!p->frozen()) {
                            // Exclude frozen particles.
                            auto r = p->position();
                            r = this->bc_->apply(r);
                            p->position(r);
                        }
                    }
                });
            }

            if ( counter % nWrite == 0) {
                // Log some simulation data.
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
                if (counter % removeCenterOfMassMotionAfter == 0) {
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
        interactor_->complete();
        logger.info("Simulation completed.");

        logger.trace("Leaving.");
    }

}

