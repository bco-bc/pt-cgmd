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
 * File:   sim-data.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:54 PM
 */

#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP

#include "stypes.hpp"
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
        energy_t ekin;
        
        /**
         * Potential energy.
         */
        energy_t epot;
        
        /**
         * Temperature.
         */
        temperature_t temperature;
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

