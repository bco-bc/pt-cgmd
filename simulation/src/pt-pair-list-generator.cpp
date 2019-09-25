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
 * File:   pt-pair-list-generator.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 25, 2019, 1:25 PM
 */

#include "simploce/simulation/pt-pair-list-generator.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/protonatable-bead.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/util/cvector_t.hpp"

namespace simploce {
    
    ProtonTransferPairListGenerator::ProtonTransferPairListGenerator(const length_t& rmax,
                                                                     const bc_ptr_t& bc) :
        rmax_{rmax}, bc_{bc}
    {        
    }
    
    ProtonTransferPairListGenerator::prot_bead_pair_list_t 
    ProtonTransferPairListGenerator::generate(const cg_ptr_t& cg) const
    {
        real_t rmax2 = rmax_() * rmax_();
        bc_ptr_t bc = bc_;
        
        return cg->doWithProtBeads<prot_bead_pair_list_t>([bc, rmax2] (const std::vector<prot_bead_ptr_t>& beads) {
            prot_bead_pair_list_t pairlist{};
            for (std::size_t i = 0; i != beads.size() - 1; ++i) {
                prot_bead_ptr_t beadi = beads[i];
                const position_t& ri = beadi->position();
                for (std::size_t j = i + 1; j != beads.size(); ++j) {
                    prot_bead_ptr_t beadj = beads[j];
                    const position_t& rj = beadj->position();
                    dist_vect_t R = bc->apply(ri, rj);
                    real_t R2 = norm2<real_t>(R);
                    if ( R2 < rmax2 ) {
                        prot_bead_pair_t pair = std::make_pair(beadi, beadj);
                        pairlist.push_back(pair);
                    }
                }
            }
            return pairlist;
        });        
    }
}
