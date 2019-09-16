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
 * File:   pol-water-force-field.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 30, 2019, 2:17 PM
 */

#include "simploce/simulation/cg-pol-water.hpp"
#include "simploce/simulation/sconf.hpp"

namespace simploce {
    
    // Interaction parameters, from Riniker et al, 2011.
    static const real_t EPS_R = 78.5;             // Relative permittivity
    static const length_t R_CW_DP = 0.2;          // nm.
    static particle_spec_ptr_t CW{};
    static particle_spec_ptr_t DP{};
    static const real_t FC = 2.0e+06;             // Force constant in kJ/(mol nm^4)
    static const real_t C12_CW_CW = 1.298e-03;    // kJ/(mol nm^12)
    static const real_t C6_CW_CW = 0.088;         // kJ/(mol nm^6)
    
    CoarseGrainedPolarizableWater::CoarseGrainedPolarizableWater(const bc_ptr_t& bc) :
        CoarseGrainedForceField{}, bc_{bc}
    {        
    }
    
    energy_t CoarseGrainedPolarizableWater::interact(const std::vector<bead_ptr_t>& all,
                                                     const std::vector<bead_ptr_t>& free,
                                                     const std::vector<bead_group_ptr_t>& groups,
                                                     const std::vector<bead_pair_list_t>& pairLists)
    {
        return 0.0;
    }
    
    std::string CoarseGrainedPolarizableWater::id() const
    {
        return conf::POLARIZABLE_WATER;
    }
    
    length_t CoarseGrainedPolarizableWater::idealDistanceCWDP()
    {
        return R_CW_DP;
    }
}
