/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 20, 2019, 12:39 PM
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "simploce/types/s-types.hpp"
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
         * @param param Simulation parameters.
         * @param particleSystem Particular system.
         * @param displacer displacer.
         */
        Simulation(param_ptr_t param,
                   p_system_ptr_t particleSystem,
                   displacer_ptr_t displacer);
        
        /**
         * Performs the simulation.
         * @param particles Output trajectory stream.
         * @param dataStream Output simulation data stream for monitoring the simulation.
         */
        void perform(std::ofstream& particles,
                     std::ofstream& dataStream);
        
    private:
        
        param_ptr_t param_;
        p_system_ptr_t particleSystem_;
        displacer_ptr_t displacer_;

    };


    
}

#endif /* SIMULATION_HPP */

