/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:56 PM
 */

#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/s-conf.hpp"
#include <iomanip>
#include <iostream>

namespace simploce {
    
    SimulationData::SimulationData() :
        t{0.0}, kinetic{0.0}, bonded{0.0}, nonBonded{0.0}, temperature{0.0},
        pressure{0.0}, numberOfProtonTransferPairs{0}, accepted{false}, acceptanceRatio{0.0} {
    }
        
    std::ostream& operator << (std::ostream& stream, const SimulationData& data) {
        const int precision = conf::PRECISION;
        const char space = conf::SPACE;
        
        stream.setf(std::ios::scientific);
        stream.precision(precision);
        energy_t total = data.kinetic + data.bonded + data.nonBonded;
        stream << std::setw(conf::REAL_WIDTH) << data.t
               << space << std::setw(conf::REAL_WIDTH) << data.kinetic
               << space << std::setw(conf::REAL_WIDTH) << data.bonded
               << space << std::setw(conf::REAL_WIDTH) << data.nonBonded
               << space << std::setw(conf::REAL_WIDTH) << total
               << space << std::setw(conf::REAL_WIDTH) << data.temperature
               << space << std::setw(conf::REAL_WIDTH) << data.pressure
               << space << std::setw(conf::REAL_WIDTH) << data.numberOfProtonTransferPairs
               << space << data.accepted
               << space << std::setw(conf::REAL_WIDTH) << data.acceptanceRatio;
#ifdef _DEBUG
        if ( total() > conf::LARGE ) {
            std::clog << "Total energy: "  << total << std::endl;
            throw std::domain_error("Very high total energy.");
        }
#endif
        return stream;
    }
    
}

