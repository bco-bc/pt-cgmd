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
 * File:   pt-pair-list-generator.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 25, 2019, 1:05 PM
 */

#ifndef PT_PAIR_LIST_GENERATOR_HPP
#define PT_PAIR_LIST_GENERATOR_HPP

#include "stypes.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    /**
     * Finds pairs of beads that may be involved in transferring protons. This 
     * applies to protonatables undergoing continuous changes in protonation state.
     */
    class ProtonTransferPairListGenerator {
    public:
                
        /**
         * Protonatable pair type.
         */
        using prot_pair_t = std::pair<cprot_bead_ptr_t, cprot_bead_ptr_t>;
        
        /**
         * Protonatable bead pair list type.
         */
        using prot_pair_list_t = std::vector<prot_pair_t>;
        
        /**
         * Constructor 
         * @param rmax Maximum distance for two protonatable beads to be involved
         * in proton transfer.
         * @param bc Boundary condition.
         */
        ProtonTransferPairListGenerator(const length_t& rmax,
                                        const bc_ptr_t& bc);
        
        /**
         * Generates protonatable bead pair list.
         * @param cg Coarse grained particle model.
         * @return List of pairs of beads possibly involved in proton transfer.
         */
        prot_pair_list_t generate(const cg_ptr_t& cg) const;
        
    private:
        
        length_t rmax_;
        bc_ptr_t bc_;
        
    };
    
     
}

#endif /* PT_PAIR_LIST_GENERATOR_HPP */

