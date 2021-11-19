/*
 * File:   displacer.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 13, 2019, 3:30 PM
 */

#ifndef DISPLACER_HPP
#define DISPLACER_HPP

#include "s-types.hpp"
#include "sim-data.hpp"

namespace simploce {
    
    /**
     * Displaces or changes the state of simulation model.
     */
    struct Displacer {
        
        virtual ~Displacer() = default;

        /**
         * Displaces particle system.
         * @param simulationData Simulation parameters.
         * @param particleSystem
         * @return Simulation data.
         */
        virtual SimulationData displace(const p_system_ptr_t& particleSystem) const = 0;
        
    };
}

#endif /* DISPLACER_HPP */

