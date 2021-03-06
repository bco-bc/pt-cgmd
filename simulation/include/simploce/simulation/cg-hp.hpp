/*
 * The MIT License
 *
 * Copyright 2019 ajuffer.
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
 * File:   cg-hp.hpp
 * Author: ajuffer
 *
 * Created on November 6, 2019, 1:13 PM
 */

#ifndef CG_HP_HPP
#define CG_HP_HPP

#include "cg-forcefield.hpp"

namespace simploce {
    
    /**
     * Applies to all group's bonds an harmonic potential. No other interactions.
     */
    class HarmonicPotential : public CoarseGrainedForceField {
    public:
        
        /**
         * Constructor.
         * @param catalog Particle catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @param fc Force constant harmonic potential.
         * @param Rref Reference particle distance,
         */
        HarmonicPotential(const spec_catalog_ptr_t& catalog,
                          const bc_ptr_t& bc,
                          const box_ptr_t& box,
                          real_t fc,
                          const length_t& Rref);
        
        std::pair<energy_t, energy_t> 
        interact(const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups,
                 const PairLists<Bead>& pairLists) override;
        
        energy_t bonded(const std::vector<bead_ptr_t>& all,
                        const std::vector<bead_ptr_t>& free,
                        const std::vector<bead_group_ptr_t>& groups,
                        const PairLists<Bead>& pairLists) override;
        
        std::pair<energy_t, energy_t> 
        interact(const bead_ptr_t& bead,
                 const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups) override;
        
        std::string id() const override;
        
        std::pair<lj_params_t, el_params_t> parameters() const override;
        
        
    private:
        
        spec_catalog_ptr_t catalog_;
        bc_ptr_t bc_;
        box_ptr_t box_;
        real_t fc_;
        length_t Rref_;
        
    };
}

#endif /* CG_HP_HPP */

