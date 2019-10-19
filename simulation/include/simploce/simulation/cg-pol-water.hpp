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
 * File:   pol-water-force-field.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 30, 2019, 2:02 PM
 */

#ifndef CG_POL_WATER_HPP
#define CG_POL_WATER_HPP

#include "cg-forcefield.hpp"
#include "sim-data.hpp"
#include "stypes.hpp"
#include <vector>

namespace simploce {
    
  /**
   * Force field for the polarizable coarse grained model according to 
   * Riniker and van Gunsteren, J. Chem. Phys. 134, 084119, 2011.
   */
    class CoarseGrainedPolarizableWater : public CoarseGrainedForceField {
    public:
        
        CoarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog,
                                      const bc_ptr_t& bc,
                                      bool protonatable);
        
        energy_t interact(const std::vector<bead_ptr_t>& all,
                          const std::vector<bead_ptr_t>& free,
                          const std::vector<bead_group_ptr_t>& groups,
                          const std::vector<bead_pair_list_t>& pairLists) override;
        
        energy_t bonded(const std::vector<bead_ptr_t>& all,
                        const std::vector<bead_ptr_t>& free,
                        const std::vector<bead_group_ptr_t>& groups,
                        const std::vector<bead_pair_list_t>& pairLists) override; 
        
        energy_t interact(const bead_ptr_t& bead,
                          const std::vector<bead_ptr_t>& all,
                          const std::vector<bead_ptr_t>& free,
                          const std::vector<bead_group_ptr_t>& groups) override;
        
        std::string id() const override;
        
        /**
         * Includes C12, and C6 LJ parameters for the particle specifications 
         * pair (CW, CW).
         * @return Parameters.
         */
        std::pair<lj_params_t, el_params_t> parameters() const override;
        
        /**
         * Ideal distance between CW and DP.
         * @return Distance, in nm.
         */
        static length_t idealDistanceCWDP();
        
    private:
        
        spec_catalog_ptr_t catalog_;
        bc_ptr_t bc_;
        
    };
}

#endif /* POL_WATER_FORCE_FIELD_HPP */

