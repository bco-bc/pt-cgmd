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
 * File:   interactor.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 12:31 PM
 */

#ifndef INTERACTOR_HPP
#define INTERACTOR_HPP

#include "sim-data.hpp"
#include <memory>
#include <vector>
#include <utility>
#include <string>

namespace simploce {
    
    /**
     * "One that interacts".
     * @param P Particle type.
     */
    template <typename P>
    class Interactor;
    
    /*
     * Specialization for atomistic particle model.
     */
    template <>
    class Interactor<Atom> {
    public:
        
        /**
         * Constructor
         * @param forcefield Coarse grained force field.
         * @param box SImulation box.
         * @param bc Boundary condition.
         */
        Interactor(const at_ff_ptr_t& forcefield,
                   const at_ppair_list_gen_ptr_t& pairListGenerator);
        
        /**
         * Computes force on atoms.
         * @param param Simulation parameters. Must hold "npairlists", number of steps
         * between updating particle pair lists.
         * @param at Atomistic particle model.
         * @return Potential energy.
         */
        energy_t interact(const sim_param_t& param, const at_ptr_t& at);
        
        /**
         * Calculates interaction energy of given atom with all other particles.
         * @param atom Atom.
         * @param param Simulation parameters.
         * @param at Atomistic particle model.
         * @return Interaction (potential) energy.
         */
        energy_t interact(const atom_ptr_t& atom,
                          const sim_param_t& param, 
                          const at_ptr_t& at);
        
        /**
         * Returns identifying name.
         * @return Identifying name.
         */
        std::string id() const;
                
    private:
        
        void updatePairlists_(const at_ptr_t& at);
        
        at_ff_ptr_t forcefield_;
        at_ppair_list_gen_ptr_t pairListGenerator_;
    };
    
    /**
     * Specialization for coarse grained model.
     */
    template <>
    class Interactor<Bead> {
    public:
        
        /**
         * Constructor
         * @param forcefield Coarse grained force field.
         * @param box SImulation box.
         * @param bc Boundary condition.
         */
        Interactor(const cg_ff_ptr_t& forcefield,
                   const cg_ppair_list_gen_ptr_t& pairListGenerator);
        
        /**
         * Computes force on beads.
         * @param param Simulation parameters.
         * @param cg Coarse grained particle model.
         * @return Potential energy.
         */
        energy_t interact(const sim_param_t& param, const cg_ptr_t& cg);
        
        /**
         * Calculates interaction energy of given bead with all other beads.
         * @param bead Bead.
         * @param param Simulation parameters.
         * @param at Coarse grained particle model.
         * @return Interaction (potential) energy.
         */
        energy_t interact(const bead_ptr_t& bead,
                          const sim_param_t& param, 
                          const cg_ptr_t& cg);
        /**
         * Returns identifying name.
         * @return Identifying name.
         */
        std::string id() const;
        
    private:
        
        void updatePairlists_(const cg_ptr_t& cg);
        
        cg_ff_ptr_t forcefield_;
        cg_ppair_list_gen_ptr_t pairListGenerator_;
    };
    
}

#endif /* INTERACTOR_HPP */

