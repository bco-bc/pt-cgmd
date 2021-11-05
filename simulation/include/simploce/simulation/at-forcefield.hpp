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
 * File:   at-forcefield.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 4:42 PM
 */

#ifndef AT_FORCEFIELD_HPP
#define AT_FORCEFIELD_HPP

#include "forcefield.hpp"
#include "simploce/particle/p-types.hpp"
#include <vector>

namespace simploce {
    
    /**
     * Interface for coarse grained force fields.
     */
    struct AtomisticForceField : public ForceField {
        
        /**
         * Atom pair type.
         */
        using bead_pair_t = std::pair<atom_ptr_t, atom_ptr_t>;
        
        /**
         * Bead pair lists type.
         */
        using atom_pair_list_t = std::vector<bead_pair_t>;

        virtual ~AtomisticForceField() {}
        
        /**
         * Computes forces on beads.
         * @param all All atoms.
         * @param free Free atoms.
         * @param groups All atoms groups.
         * @param pairLists One or multiple atom pair lists.
         * @return Bonded and non-bonded potential energy.
         */
        virtual std::pair<energy_t, energy_t>
        interact(const std::vector<atom_ptr_t>& all,
                 const std::vector<atom_ptr_t>& free,
                 const std::vector<atom_group_ptr_t>& groups,
                 const std::vector<atom_pair_list_t>& pairLists) = 0;
    };
}


#endif /* AT_FORCEFIELD_HPP */

