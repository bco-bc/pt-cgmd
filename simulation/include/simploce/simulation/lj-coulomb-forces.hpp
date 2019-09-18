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
 * File:   lj-coulomb.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 17, 2019, 3:06 PM
 */

#ifndef LJ_COULOMB_HPP
#define LJ_COULOMB_HPP

#include "cg-forcefield.hpp"
#include "stypes.hpp"
#include "simploce/util/map2.hpp"
#include <string>
#include <map>
#include <vector>

namespace simploce {
    
    /**
     * Calculates LJ and Coulomb interaction.
     * @param P Particle type.
     */
    template <typename P>
    class LJCoulombForces;
    
    /**
     * Specialization for beads.
     */
    template <>
    class LJCoulombForces<Bead> : public CoarseGrainedForceField {
    public:
        
        /**
         * Type for holding pairs (C12, C6) of LJ interaction parameters 
         * for pairs of particle specs, each identified by a specification name. 
         */
        using lj_params_t = MatrixMap<std::string, std::pair<real_t, real_t>>;
        
        /**
         * Type for holding parameters for electrostatic interactions.
         */
        using el_params_t = std::map<std::string, real_t>;
        
        LJCoulombForces(const lj_params_t& ljParams, 
                        const el_params_t& elParams,
                        const bc_ptr_t& bc);
        
        energy_t interact(const std::vector<bead_ptr_t>& all,
                          const std::vector<bead_ptr_t>& free,
                          const std::vector<bead_group_ptr_t>& groups,
                          const std::vector<bead_pair_list_t>& pairLists) override;
        
        std::string id() const override;
        
    private:
        
        lj_params_t ljParams_;
        el_params_t elParams_;
        bc_ptr_t bc_;
    };
}

#endif /* LJ_COULOMB_HPP */

