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
 * File:   pt-langevin-velocity-verlet.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 25, 2019, 2:41 PM
 */

#ifndef PT_LANGEVIN_VELOCITY_VERLET_HPP
#define PT_LANGEVIN_VELOCITY_VERLET_HPP

#include "cg-displacer.hpp"
#include "s-types.hpp"

namespace simploce {
    
    /**
     * Displaces protonatable beads with continuous varying charges and mass values.
     */
    class ProtonTransferLangevinVelocityVerlet : public CoarseGrainedDisplacer {
    public:
        
        ProtonTransferLangevinVelocityVerlet(const cg_interactor_ptr_t& interactor,
                                             const prot_pair_list_gen_ptr_t& generator,
                                             const pt_displacer_ptr_t& displacer);
        
        SimulationData displace(const sim_param_t& param, 
                                const cg_mod_ptr_t& cg) const override;
        
        std::string id() const override;
        
    private:
        
        cg_interactor_ptr_t interactor_;
        prot_pair_list_gen_ptr_t generator_;
        pt_displacer_ptr_t displacer_;
    };
}

#endif /* PT_LANGEVIN_VELOCITY_VERLET_HPP */

