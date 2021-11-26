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

    Simulation::Simulation(sim_param_ptr_t simulationParams,
                           p_system_ptr_t particleSystem,
                           displacer_ptr_t displacer) :
        simulationParameters_{std::move(simulationParams)},
        particleSystem_{std::move(particleSystem)},
        displacer_{std::move(displacer)} {
    }

    void Simulation::perform(std::ofstream& trajectoryStream,
                             std::ofstream& dataStream) {
        static util::Logger logger("simploce::simulation::perform()");

        std::size_t numberAccepted = 0;
        auto box = particleSystem_->box();
        auto nSteps = simulationParameters_->get<int>("simulation.nsteps", 10000);
        auto nWrite = simulationParameters_->get<int>("simulation.nwrite", 10);
        for (int counter = 1; counter <= nSteps; ++counter) {
            logger.trace("Step #: " + util::toString(counter));
            SimulationData data = displacer_->displace(particleSystem_);
            if ( data.accepted ) {
                numberAccepted += 1;
            }
            if ( counter % nWrite == 0) {
                data.pressure =
                        particleSystem_->doWithAllFreeGroups<pressure_t>([box, data] (
                                const std::vector<p_ptr_t>& all,
                                const std::vector<p_ptr_t>& free,
                                const std::vector<pg_ptr_t>& groups) {
                    return properties::pressure(all, data.temperature, box);
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

