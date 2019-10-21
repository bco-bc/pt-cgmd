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
 * File:   cg-forcefield.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 4:13 PM
 */

#ifndef CG_FORCEFIELD_HPP
#define CG_FORCEFIELD_HPP

#include "forcefield.hpp"
#include "pair-lists.hpp"
#include "stypes.hpp"
#include <vector>

namespace simploce {
    
    /**
     * Interface for coarse grained force fields.
     */
    struct CoarseGrainedForceField : public ForceField {
        
        /**
         * Bead pair type.
         */
        using bead_pair_t = std::pair<bead_ptr_t, bead_ptr_t>;
        
        /**
         * Bead pair lists type.
         */
        using bead_pair_list_t = std::vector<bead_pair_t>;
                
        virtual ~CoarseGrainedForceField() {}
        
        /**
         * Computes forces due to -all- interactions on beads. Updates/adds all forces 
         * acting on beads.
         * @param all All beads.
         * @param free Free beads.
         * @param groups All bead groups.
         * @param pairLists One or multiple bead pair lists.
         * @return Potential energy.
         */
        virtual energy_t interact(const std::vector<bead_ptr_t>& all,
                                  const std::vector<bead_ptr_t>& free,
                                  const std::vector<bead_group_ptr_t>& groups,
                                  const PairLists<Bead>& pairLists) = 0;
        
        /**
         * Computes forces due to bonded interactions on beads. Updates/adds forces 
         * bonded acting on beads.
         * @param all All beads.
         * @param free Free beads.
         * @param groups  All bead groups.
         * @param pairLists  One or multiple bead pair lists.
         * @return Potential energy.
         */
        virtual energy_t bonded(const std::vector<bead_ptr_t>& all,
                                const std::vector<bead_ptr_t>& free,
                                const std::vector<bead_group_ptr_t>& groups,
                                const PairLists<Bead>& pairLists) = 0;
        
        
        /**
         * Returns interaction energy of given bead with all other bead.
         * @param bead Bead
         * @param all All beads.
         * @param free Free beads.
         * @param groups bead groups.
         * @return Energy.
         */
        virtual energy_t interact(const bead_ptr_t& bead,
                                  const std::vector<bead_ptr_t>& all,
                                  const std::vector<bead_ptr_t>& free,
                                  const std::vector<bead_group_ptr_t>& groups) = 0;
        
    };
}


#endif /* CG_FORCEFIELD_HPP */

