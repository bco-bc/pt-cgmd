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
#include "simploce/particle/continuous-protonatable-bead.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/util/cvector_t.hpp"

namespace simploce {
    
    ProtonTransferPairListGenerator::ProtonTransferPairListGenerator(const length_t& rmax,
                                                                     const bc_ptr_t& bc) :
        rmax_{rmax}, bc_{bc}
    {        
    }
    
    ProtonTransferPairListGenerator::prot_pair_list_t 
    ProtonTransferPairListGenerator::generate(const cg_ptr_t& cg) const
    {
        real_t rmax2 = rmax_() * rmax_();
        bc_ptr_t bc = bc_;
        
        return cg->doWithProtBeads<prot_pair_list_t>([bc, rmax2] (const std::vector<dprot_bead_ptr_t>& discrete,
                                                                  const std::vector<cprot_bead_ptr_t>& continuous) {
            prot_pair_list_t pairlist{};
            if ( !continuous.empty() ) {
                for (auto i = continuous.begin(); i != (continuous.end() - 1); ++i) {
                    auto beadi = *i;
                    const auto& ri = beadi->position();
                    for (auto j = i + 1; j != continuous.end(); ++j) {
                        auto beadj = *j;
                        const auto& rj = beadj->position();
                        auto R = bc->apply(ri, rj);
                        real_t R2 = norm2<real_t>(R);
                        if ( R2 < rmax2 ) {
                            prot_pair_t pair = std::make_pair(beadi, beadj);
                            pairlist.push_back(pair);
                        }
                    }
                }
            }
            return pairlist;
        });        
    }
}
