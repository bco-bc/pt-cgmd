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
 * File:   simulation.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 10, 2019, 1:53 PM
 */

#include "simploce/simulation/simulation.hpp"
#include "simploce/simulation/sim-model.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/sconf.hpp"
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace simploce {
    
    Simulation<Bead>::Simulation(const cg_sim_model_ptr_t& sm) : sm_{sm}
    {        
    }
    
    void Simulation<Bead>::perform(const sim_param_t& param,
                                   std::ofstream& trajStream,
                                   std::ofstream& dataStream)
    {
        const auto width = conf::WIDTH;
        const auto space = conf::SPACE;
        
        if ( sm_->size() == 0 ) {
            throw std::domain_error(
                "No particles! Nothing to simulate."
            );
        }

        std::size_t nsteps = param.get<std::size_t>("nsteps", 10000);
        std::size_t nwrite = param.get<std::size_t>("nwrite", 10);
        
        for (std::size_t counter = 1; counter <= nsteps; ++counter) {
            SimulationData data = sm_->displace(param);
            if ( counter % nwrite == 0 ) {
                dataStream << std::setw(width) << counter << space << data << std::endl;
                dataStream.flush();
                sm_->saveState(trajStream);
            }
        }
                
    }
    
}
