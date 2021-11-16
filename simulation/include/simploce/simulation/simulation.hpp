/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 20, 2019, 12:39 PM
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "s-types.hpp"
#include "interactor.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <memory>
#include <fstream>
#include <iostream>

namespace simploce {
    
    /**
     * Performs a simulation. Writes particle state and simulation data to separate
     * output streams.
     * @tparam p Particle type.
     */
    template <typename P>
    class Simulation {
    public:

        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;

        /**
         * Particle system pointer type.
         */
        using p_sys_ptr_t = std::shared_ptr<ParticleSystem<P>>;

        /**
         * Interactor pointer type.
         */
        using interactor_ptr_t = std::shared_ptr<Interactor<P>>;

        /**
         * Constructor. All arguments are required.
         * @param Simulation parameters.
         */
        Simulation(sim_param_ptr_t simulationParams,
                   particle_system_fact_ptr_t particleSystem,
                   interactor_ptr_t interactor);
        
        /**
         * Performs the simulation.
         * @param trajStream Output trajectory stream.
         * @param dataStream Output simulation data stream for monitoring the simulation.
         */
        void perform(std::ofstream& trajStream,
                     std::ofstream& dataStream);
        
    private:
        
        sim_param_ptr_t simulationParameters_;
        particle_system_fact_ptr_t particleSystem_;
        interactor_ptr_t interactor_;

    };

    template <typename P>
    Simulation<P>::Simulation(sim_param_ptr_t simulationParams,
                              particle_system_fact_ptr_t particleSystem,
                              interactor_ptr_t interactor) :
        simulationParameters_{std::move(simulationParams)},
        particleSystem_{std::move(particleSystem)},
        interactor_{std::move(interactor)} {
    }

    template <typename P>
    void Simulation<P>::perform(std::ofstream& trajStream,
                                std::ofstream& dataStream) {
        static util::Logger logger("simploce::simulation::perform()");

        auto box = particleSystem->box();
        auto nSteps = param.get<int>("numberOfSteps", 10000);
        auto nWrite = param.get<int>("nwrite", 10);
        for ( counter = 1; counter != nSteps; ++counter) {
            logger.debug("Step #: " + util::toString(counter));
            SimulationData data; // From displacer...
            if ( counter % nWrite == 0) {
                data.pressure =
                        particleSystem_->doWithAllFreeGroups<pressure_t>([box, data] (
                                const std::vector<p_ptr_t>& all,
                                const std::vector<p_ptr_t>& free,
                                const std::vector<p_group_ptr_t>& groups) {
                    return properties::pressure(all, data.temperature, box);
                });
                dataStream << std::setw(conf::INTEGER_WIDTH) << counter << space << data << std::endl;
                particleSystem_->saveState(trajStream);
                trajStream.flush();
                dataStream.flush();
            }
        }
    }

    
}

#endif /* SIMULATION_HPP */

