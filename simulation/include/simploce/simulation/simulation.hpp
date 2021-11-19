/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 20, 2019, 12:39 PM
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "s-types.hpp"
#include <fstream>

namespace simploce {
    
    /**
     * Performs a simulation. Writes particle state and simulation data to separate
     * output streams.
     */
    class Simulation {
    public:

        /**
         * Constructor. All arguments are required.
         * @param simulationParams Simulation parameters.
         * @param particleSystem Particular system.
         * @param displacer Displacer.
         */
        Simulation(sim_param_ptr_t simulationParams,
                   p_system_ptr_t particleSystem,
                   displacer_ptr_t displacer);
        
        /**
         * Performs the simulation.
         * @param trajectoryStream Output trajectory stream.
         * @param dataStream Output simulation data stream for monitoring the simulation.
         */
        void perform(std::ofstream& trajectoryStream,
                     std::ofstream& dataStream);
        
    private:
        
        sim_param_ptr_t simulationParameters_;
        p_system_ptr_t particleSystem_;
        displacer_ptr_t displacer_;

    };


    
}

#endif /* SIMULATION_HPP */

