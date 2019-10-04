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
 * File:   constant-rate-pt.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 4, 2019, 2:16 PM
 */

#include "simploce/simulation/constant-rate-pt.hpp"
#include "simploce/particle/continuous-protonatable-bead.hpp"
#include <stdexcept>

namespace simploce {
    
    using prot_pair_t = ProtonTransferPairListGenerator::prot_pair_t;
    
    static bool hasProton_(const prot_pair_t& pair)
    {
        return pair.first->isProtonated() || pair.second->isProtonated();
    }
    
    ConstantRateProtonTransfer::ConstantRateProtonTransfer(const rate_t& rate) :
        rate_{rate}
    {   
        if ( rate_() < 0.0 ) {
            throw std::domain_error("A rate value must be nonnegative number.");
        }
    }
        
    void 
    ConstantRateProtonTransfer::transfer(const sim_param_t& param,
                                         const prot_pair_list_t& pairList) const
    {
        for (auto pair : pairList) {
            if ( hasProton_(pair) ) {
                auto p1 = pair.first;
                auto ps1 = p1->protonationState();                
                auto p2 = pair.second;
                auto ps2 = p2->protonationState();
                
                // Transfer only if one of the beads is protonated.
                if (ps1 == 1 && ps2 == 0 ) {
                    
                } else if ( ps1 == 0 && ps2 == 1) {
                    
                }
            }
        }
    }
    
}