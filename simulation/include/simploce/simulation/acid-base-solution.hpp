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
 * File:   acid-base-solution.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 8, 2019, 3:28 PM
 */

#ifndef ACID_BASE_SOLUTION_HPP
#define ACID_BASE_SOLUTION_HPP

#include "cg-forcefield.hpp"
#include "sim-data.hpp"
#include "stypes.hpp"
#include <vector>

namespace simploce {
    
    /**
     * Coarse grained protonatable acids and/or bases in coarse grained polarizable water.
     */
    class AcidBaseSolution : public CoarseGrainedForceField {
    public:
        
        /**
         * Constructor.
         * @param catalog Particle specifications catalog.
         * @param bc Boundary conditions.
         * @param box SImulation box.
         * @param water Protonatable water.
         */
        AcidBaseSolution(const spec_catalog_ptr_t& catalog,
                         const bc_ptr_t& bc,
                         const box_ptr_t& box,
                         const cg_ff_ptr_t& water);
        
        energy_t interact(const std::vector<bead_ptr_t>& all,
                          const std::vector<bead_ptr_t>& free,
                          const std::vector<bead_group_ptr_t>& groups,
                          const PairLists<Bead>& pairLists) override;
        
        energy_t bonded(const std::vector<bead_ptr_t>& all,
                        const std::vector<bead_ptr_t>& free,
                        const std::vector<bead_group_ptr_t>& groups,
                        const PairLists<Bead>& pairLists) override;
        
        energy_t interact(const bead_ptr_t& bead,
                          const std::vector<bead_ptr_t>& all,
                          const std::vector<bead_ptr_t>& free,
                          const std::vector<bead_group_ptr_t>& groups) override;
        
        std::string id() const override;
        
        std::pair<lj_params_t, el_params_t> parameters() const override;
        
    private:
        
        spec_catalog_ptr_t catalog_;
        bc_ptr_t bc_;
        box_ptr_t box_;
        cg_ff_ptr_t water_;
        
    };
}

#endif /* ACID_BASE_SOLUTION_HPP */

