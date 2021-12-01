/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:54 PM
 */

#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP

#include "s-types.hpp"
#include <iostream>

namespace simploce {
    
    struct SimulationData {
        
        SimulationData();
        
        /**
         * Time
         */
        stime_t t;
        
        /**
         * Kinetic energy.
         */
        energy_t kinetic;
        
        /**
         * Bonded potential energy.
         */
        energy_t bonded;
        
        /**
         * Non bonded potential energy.
         */
        energy_t nonBonded;

        /**
         * Energy due to external potentials.
         */
        energy_t external;

        /**
         * Total momentum (norm of the vector).
         */
        real_t totalMomentum;
        
        /**
         * Temperature.
         */
        temperature_t temperature;
        
        /**
         * Pressure.
         */
        pressure_t pressure;
        
        /**
         * Number of particle pairs possibly involved in proton transfer.
         */
        std::size_t numberOfProtonTransferPairs;
        
        /**
         * Move was accepted in a Monte Carlo simulation.
         */
        bool accepted;
        
        /**
         * Acceptance ratio (%) in a Monte Carlo simulation.
         */
        real_t acceptanceRatio;
    };
    
    /**
     * Writes simulation data to output stream.
     * @param stream Output stream.
     * @param data Simulation data.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const SimulationData& data);
}

#endif /* SIM_DATA_HPP */

