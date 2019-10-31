/*
 * The MIT License
 *
 * Copyright 2019 juffer.
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

#include "simploce/analysis/analysis.hpp"
#include "simploce/analysis/analyzer.hpp"
#include "simploce/simulation/sim-model.hpp"
#include <iostream>

namespace simploce {
    
    Analysis<Bead>::Analysis(const cg_sim_model_ptr_t& sm,
                             const cg_analyzer_ptr_t& analyzer) :
        sm_{sm}, analyzer_{analyzer}
    {        
    }
        
    void 
    Analysis<Bead>::perform(const sim_param_t& param,
                            std::istream& trajectory)
    {
        const auto nskip = param.get<std::size_t>("nskip");
        std::clog << "Skipping first " << nskip << " states in trajectory." << std::endl;
        
        std::size_t counter = 0;
        sm_->readState(trajectory);
        while (trajectory.good() ) {
            counter += 1;
            if ( counter > nskip ) {
                sm_->doWithAllFreeGroups<void>([this] (const std::vector<bead_ptr_t>& all,
                                                       const std::vector<bead_ptr_t>& free,
                                                       const std::vector<bead_group_ptr_t>& groups) {
                    this->analyzer_->perform(all, free, groups);
                });
            }
            sm_->readState(trajectory);
        }
    }
    
}