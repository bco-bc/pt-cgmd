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
 * File:   interactor.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 12:49 PM
 */

#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/simulation/at-forcefield.hpp"
#include "simploce/simulation/cg-forcefield.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    // Atom pair type
    using atom_pair_t = std::pair<atom_ptr_t, atom_ptr_t>;
        
    //Bead pair lists type (single list).
    using atom_pair_list_t = std::vector<atom_pair_t>;
        
    
    // Bead pair type
    using bead_pair_t = std::pair<bead_ptr_t, bead_ptr_t>;
        
    //Bead pair lists type (single list).
    using bead_pair_list_t = std::vector<bead_pair_t>;
    
    // Bead pair lists.
    static std::vector<atom_pair_list_t> atomPairLists_{};

    // Bead pair lists.
    static std::vector<bead_pair_list_t> beadPairLists_{};
    
    
    Interactor<Atom>::Interactor(const at_ff_ptr_t& forcefield,
                                 const at_ppair_list_gen_ptr_t& pairListGenerator) :
        forcefield_{forcefield}, pairListGenerator_{pairListGenerator}
    {        
    }
        
    energy_t 
    Interactor<Atom>::interact(const sim_param_t& param, 
                               const at_ptr_t& at)
    {
        static bool setup = false;
        static std::size_t npairlist = 0;
        static std::size_t counter = 0;
        
        if ( !setup ) {
            npairlist = param.get<std::size_t>("npairlists");
            setup = true;
        }
        if ( counter % npairlist == 0) {
            this->updatePairlists_(at);
        }
        
        energy_t epot = at->doWithAllFreeGroups<energy_t>([this] (const std::vector<atom_ptr_t>& all,
                                                                  const std::vector<atom_ptr_t>& free,
                                                                  const std::vector<atom_group_ptr_t>& groups) {
            for (auto p: all) {
                p->resetForce();
            }
            return this->forcefield_->interact(all, free, groups, atomPairLists_);
        });
        
        counter += 1;
        return epot;
    }
    
    void 
    Interactor<Atom>::updatePairlists_(const at_ptr_t& at)
    {
        atomPairLists_ = 
            at->doWithAllFreeGroups<std::vector<atom_pair_list_t>>([this] (const std::vector<atom_ptr_t>& all,
                                                                           const std::vector<atom_ptr_t>& free,
                                                                           const std::vector<atom_group_ptr_t>& groups) {
                return this->pairListGenerator_->generate(all, free, groups);
            });
    }
    
    std::string 
    Interactor<Atom>::id() const
    {
        return forcefield_->id();
    }

    
    Interactor<Bead>::Interactor(const cg_ff_ptr_t& forcefield,
                                 const cg_ppair_list_gen_ptr_t& pairListGenerator) :
        forcefield_{forcefield}, pairListGenerator_{pairListGenerator}
    {      
    }
        
    energy_t 
    Interactor<Bead>::interact(const sim_param_t& param, 
                               const cg_ptr_t& cg)
    {
        static bool setup = false;
        static std::size_t npairlist = 0;
        static std::size_t counter = 0;
        
        if ( !setup ) {
            npairlist = param.get<std::size_t>("npairlists");
            setup = true;
        }
        if ( counter % npairlist == 0) {
            this->updatePairlists_(cg);
        }
        
        energy_t epot = 
            cg->doWithAllFreeGroups<energy_t>([this] (const std::vector<bead_ptr_t>& all,
                                                      const std::vector<bead_ptr_t>& free,
                                                      const std::vector<bead_group_ptr_t>& groups) {
                for (auto p: all) {
                    p->resetForce();
                }
                return this->forcefield_->interact(all, free, groups, beadPairLists_);
            });
        
        counter += 1;
        return epot;
    }
    
    void
    Interactor<Bead>::updatePairlists_(const cg_ptr_t& cg)
    {
        beadPairLists_ = 
            cg->doWithAllFreeGroups<std::vector<bead_pair_list_t>>([this] (const std::vector<bead_ptr_t>& all,
                                                                           const std::vector<bead_ptr_t>& free,
                                                                           const std::vector<bead_group_ptr_t>& groups) {
                return this->pairListGenerator_->generate(all, free, groups);
            });
    }
    
    std::string
    Interactor<Bead>::id() const
    {
        return forcefield_->id();
    }

}