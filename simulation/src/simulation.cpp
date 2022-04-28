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
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/util.hpp"

namespace simploce {

    Simulation::Simulation(param_ptr_t param,
                           p_system_ptr_t particleSystem,
                           displacer_ptr_t displacer) :
            param_{std::move(param)},
            particleSystem_{std::move(particleSystem)},
            displacer_{std::move(displacer)} {
    }

    void Simulation::perform(std::ofstream& trajectoryStream,
                             std::ofstream& dataStream) {
        static util::Logger logger("simploce::simulation::perform()");

        std::size_t numberAccepted = 0;
        auto box = particleSystem_->box();
        auto nSteps = param_->get<int>("simulation.nsteps", 10000);
        auto nWrite = param_->get<int>("simulation.nwrite", 10);
        logger.debug("Number of steps: " + util::toString(nSteps));
        logger.debug("Number of steps between writing simulation data: " + util::toString(nWrite));
        for (int counter = 1; counter <= nSteps; ++counter) {
            logger.trace("Step #: " + util::toString(counter));
            SimulationData data = displacer_->displace(particleSystem_);
            if ( data.accepted ) {
                numberAccepted += 1;
            }
            if ( counter % nWrite == 0) {
                // Pressure is computed for displaceable particles, not for all particles.
                data.pressure =
                        particleSystem_->doWithDisplaceables<pressure_t>([box, data] (
                                const std::vector<p_ptr_t>& particles) {
                    return properties::pressure(particles, data.temperature, box);
                });
                data.acceptanceRatio = real_t(numberAccepted)/real_t(counter) * 100.0;
                dataStream << std::setw(conf::INTEGER_WIDTH) << counter << conf::SPACE << data << std::endl;
                particleSystem_->writeState(trajectoryStream);
                trajectoryStream.flush();
                dataStream.flush();
            }
        }
    }

}

