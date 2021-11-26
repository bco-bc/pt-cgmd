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
 * File:   pt.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 4, 2019, 2:05 PM
 */

#ifndef PT_HPP
#define PT_HPP

#include "pt-pair-list-generator.hpp"
#include "s-types.hpp"
#include <vector>

namespace simploce {
    
    /**
     * Interface for transferring protons between two protonatable beads with 
     * continuously varying protonation states. 
     */
    struct ProtonTransfer {
        
        /**
         * Protonatable bead pair list typeName.
         */
        using prot_pair_list_t = ProtonTransferPairListGenerator::prot_pair_list_t;
        
        virtual ~ProtonTransfer() {}
        
        /**
         * Transfers proton between two beads.
         * @param param Simulation parameters.
         * @param continuous Protonatable beads with continuously varying 
         * protonation states.
         * @param pairList Pairs of protonatables possibly involved in proton transfer.
         */
        virtual void transfer(const sim_param_t& param,
                              const std::vector<prot_bead_ptr_t>& continuous,
                              const prot_pair_list_t& pairList) const = 0;
    };
}

#endif /* PT_HPP */

