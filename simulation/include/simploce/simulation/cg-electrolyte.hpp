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
 * File:   cg-electrolyte.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 18, 2019, 1:02 PM
 */

#ifndef CG_ELECTROLYTE_HPP
#define CG_ELECTROLYTE_HPP

#include "cg-forcefield.hpp"

namespace simploce {
    
    /**
     * Simple force field for a coarse grained electrolyte solution. Water is a background
     * continuum, ions interaction through a screened Coulomb and Lennard-Jones potential.
     * Interaction parameters are taken from 
     * Lenart et al, J. Chem. Phys., 126, 044509, 2007.
     */
    class CoarseGrainedElectrolyte : public CoarseGrainedForceField {
    public:
        
        CoarseGrainedElectrolyte(const spec_catalog_ptr_t& catalog,
                                 const bc_ptr_t& bc);
        
        energy_t interact(const std::vector<bead_ptr_t>& all,
                          const std::vector<bead_ptr_t>& free,
                          const std::vector<bead_group_ptr_t>& groups,
                          const std::vector<bead_pair_list_t>& pairLists) override;
        
        energy_t bonded(const std::vector<bead_ptr_t>& all,
                        const std::vector<bead_ptr_t>& free,
                        const std::vector<bead_group_ptr_t>& groups,
                        const std::vector<bead_pair_list_t>& pairLists) override;
        
        std::string id() const override;
        
        std::pair<lj_params_t, el_params_t> parameters() const override;
        
    private:
        
        spec_catalog_ptr_t catalog_;
        bc_ptr_t bc_;
        
    };
}

#endif /* CG_ELECTROLYTE_HPP */

