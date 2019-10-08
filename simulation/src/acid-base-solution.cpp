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
 * File:   acid-base-solution.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 8, 2019, 3:47 PM
 */

#include "simploce/simulation/acid-base-solution.hpp"
#include "simploce/simulation/sconf.hpp"
#include <stdexcept>

namespace simploce {
    
    using lj_params_t = ForceField::lj_params_t;
    using el_params_t = ForceField::el_params_t;
    
    static lj_params_t ljParams_;
    static el_params_t elParams_;
    
    void
    setup_(const spec_catalog_ptr_t& catalog,
           const bc_ptr_t& bc,
           const cg_ff_ptr_t& water)
    {
        // Water parameters.
        auto parameters = water->parameters();
        ljParams_ = parameters.first;
        elParams_ = parameters.second;
        
        // LJ
        
        
        
    }
    
    AcidBaseSolution::AcidBaseSolution(const spec_catalog_ptr_t& catalog,
                                       const bc_ptr_t& bc,
                                       const cg_ff_ptr_t& water) :
        catalog_{catalog}, bc_{bc}, water_{water}
    {   
        if ( !catalog_ ) {
            throw std::domain_error(
                "AcidBaseSolution: Missing particle specification catalog."
            );
        }
        if ( !bc_ ) {
            throw std::domain_error(
                "AcidBaseSolution: Missing boundary condition."
            );            
        }
        if ( !water ) {
            throw std::domain_error(
                "AcidBaseSolution: Missing protonatable water."
            );                        
        }
        
        setup_(catalog_, bc_, water_);
    }
    
    energy_t 
    AcidBaseSolution::interact(const std::vector<bead_ptr_t>& all,
                               const std::vector<bead_ptr_t>& free,
                               const std::vector<bead_group_ptr_t>& groups,
                               const std::vector<bead_pair_list_t>& pairLists)
    {
        return 0.0;
    }
    
    energy_t 
    AcidBaseSolution::bonded(const std::vector<bead_ptr_t>& all,
                             const std::vector<bead_ptr_t>& free,
                             const std::vector<bead_group_ptr_t>& groups,
                             const std::vector<bead_pair_list_t>& pairLists)
    {
        return water_->bonded(all, free, groups, pairLists);
    }
    
    std::string 
    AcidBaseSolution::id() const
    {
        return conf::ACID_BASE_SOLUTION;
    }
    
    std::pair<lj_params_t, el_params_t> 
    AcidBaseSolution::parameters() const
    {
        return water_->parameters();
    }
}
