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
 * File:   pt-langevin-velocity-verlet.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 25, 2019, 2:47 PM
 */

#include "simploce/simulation/pt-langevin-velocity-verlet.hpp"
#include "simploce/simulation/pt-pair-list-generator.hpp"
#include <stdexcept>

namespace simploce {
    
    ProtonTransferLangevinVelocityVerlet::
        ProtonTransferLangevinVelocityVerlet(const cg_interactor_ptr_t& interactor,
                                             const pt_pair_list_gen_ptr_t& generator) : 
    interactor_{interactor}, generator_{generator}
    {        
        if ( !interactor ) {
            throw std::domain_error("PT LangevinVelocityVerlet: Missing interactor.");
        }
        if ( !generator ) {
            throw std::domain_error(
                "PT LangevinVelocityVerlet: Missing protonatable bead pair list generator."
            );
        }
    }
    
    SimulationData ProtonTransferLangevinVelocityVerlet::displace(const sim_param_t& param, 
                                                                  const cg_ptr_t& cg) const
    {
        using prot_pair_list_t = 
            ProtonTransferPairListGenerator::prot_pair_list_t;
        
        static std::size_t counter = 0;
        static prot_pair_list_t pairlist{};
        
        SimulationData data{};
        counter += 1;
        
        std::size_t npairlists = param.get<std::size_t>("npairlists");
        if ( counter % npairlists == 0 ) {
            pairlist = generator_->generate(cg);
            data.numberOfProtonTransferPairs = pairlist.size();
        }
        
        // Transfer protons.
        if ( !pairlist.empty() ) {
            
        }
        
        // Update positions and velocities.
        
        return SimulationData{};
    }
}

