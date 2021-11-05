/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 20, 2019, 12:39 PM
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "s-types.hpp"
#include <fstream>
#include <iostream>

namespace simploce {
    
    /**
     * Performs a simulation. Writes particle state and simulation data to separate
     * output streams.
     * @param P Particle type.
     */
    template <typename P>
    class Simulation;
    
    template <>
    class Simulation<Bead> {
    public:
        
        /**
         * Constructor.
         * @param sm Simulation model.
         */
        Simulation(const cg_sim_model_ptr_t& sm);
        
        /**
         * Performs the simulation.
         * @param param Parameters. Must provide,
         * <ul>
         *  <li>nsteps: Number of steps.</li>
         *  <li>
         *      nwrite: Number of steps between writing simulation data and saving state
         *      in the trajectectory.
         *  </li>
         * </ul>
         * @param trajStream Output trajectory stream.
         * @param dataStream Output simulation data stream.
         */
        void perform(const sim_param_t& param,
                     std::ofstream& trajStream,
                     std::ofstream& dataStream);
        
    private:
        
        cg_sim_model_ptr_t sm_;
    };    
    
}

#endif /* SIMULATION_HPP */

