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
         * @param catalog Particle specifications catalog.
         * @param displacer displacer.
         * @param bc Boundary condition.
         */
        Simulation(param_ptr_t param,
                   p_system_ptr_t particleSystem,
                   spec_catalog_ptr_t catalog,
                   displacer_ptr_t displacer,
                   bc_ptr_t bc);
        
        /**
         * Performs the simulation.
         * @param all Output trajectory stream.
         * @param dataStream Output simulation data stream for monitoring the simulation.
         */
        void perform(std::ofstream& all,
                     std::ofstream& dataStream);
        
    private:
        
        param_ptr_t param_;
        p_system_ptr_t particleSystem_;
        spec_catalog_ptr_t catalog_;
        displacer_ptr_t displacer_;
        bc_ptr_t bc_;

    };


    
}

#endif /* SIMULATION_HPP */

