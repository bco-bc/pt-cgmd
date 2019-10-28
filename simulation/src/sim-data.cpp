/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   sim-data.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:56 PM
 */

#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/sconf.hpp"
#include <iomanip>

namespace simploce {
    
    SimulationData::SimulationData() :
        t{0.0}, ekin{0.0}, epot{0.0}, temperature{0.0}, pressure{0.0},
        numberOfProtonTransferPairs{0}, accepted{false}, acceptanceRatio{0.0}
    {            
    }
        
    std::ostream& operator << (std::ostream& stream, const SimulationData& data)
    {
        const int width = conf::WIDTH;
        const int precision = conf::PRECISION;
        const char space = conf::SPACE;
        
        stream.setf(std::ios::scientific);
        stream.precision(precision);
        energy_t etot = data.ekin + data.epot;
        stream << std::setw(width) << data.t
               << space << std::setw(width) << data.ekin
               << space << std::setw(width) << data.epot;
        stream << space << std::setw(width) << etot
               << space << std::setw(width) << data.temperature
               << space << std::setw(width) << data.pressure
               << space << std::setw(width) << data.numberOfProtonTransferPairs
               << space << data.accepted
               << space << std::setw(width) << data.acceptanceRatio;
        
        return stream;
    }
    
}

