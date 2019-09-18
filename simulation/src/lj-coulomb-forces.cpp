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
 * File:   lj-coulomb-forces.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 17, 2019, 3:43 PM
 */

#include "simploce/simulation/lj-coulomb-forces.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/util/util.hpp"
#include <future>

namespace simploce {
    
    using lj_params_t = LJCoulombForces<Bead>::lj_params_t;
    using el_params_t = LJCoulombForces<Bead>::el_params_t;
    using result_t = std::pair<energy_t, std::vector<force_t>>;
    using bead_pair_list_t = ParticlePairListGenerator<Bead>::p_pair_list_t;
    
    static result_t calculate_(const bead_pair_list_t& single,
                               std::size_t nbeads,
                               const lj_params_t& ljParams,
                               const el_params_t& elParams,
                               const bc_ptr_t& bc)
    {
        std::vector<force_t> forces(nbeads, force_t{});
        energy_t epot{0.0};
        return std::make_pair(epot, forces);
    }
    
    LJCoulombForces<Bead>::LJCoulombForces(const lj_params_t& ljParams, 
                                           const el_params_t& elParams, 
                                           const bc_ptr_t& bc) :
        ljParams_{ljParams}, elParams_{elParams}, bc_{bc}
    {        
    }
        
    energy_t LJCoulombForces<Bead>::interact(const std::vector<bead_ptr_t>& all,
                                             const std::vector<bead_ptr_t>& free,
                                             const std::vector<bead_group_ptr_t>& groups,
                                             const std::vector<bead_pair_list_t>& pairLists)
    {   
        static std::vector<result_t> results;
        static std::vector<std::future<result_t> > futures{};
        
        // Set up tasks for concurrent force and energy calculation.
        std::size_t nlists = pairLists.size();
        std::size_t nbeads = all.size();
        
        if ( pairLists.size() > 1 ) {
            std::size_t ntasks = pairLists.size() - 1;  // One other for current thread.
            futures.clear();
            results.clear();
            for (std::size_t k = 0; k != ntasks; ++k) {
                const bead_pair_list_t& single = pairLists[k];
                futures.push_back(std::async(
                    std::launch::async, 
                    calculate_,
                    std::ref(single),
                    nbeads,
                    std::ref(ljParams_),
                    std::ref(elParams_),
                    std::ref(bc_)
                ));

            }
            // Wait for these tasks to complete.
            results = util::waitForAll<result_t>(futures);
        }
    
        // Current thread.
        const bead_pair_list_t& last = pairLists[nlists - 1];
        result_t result = calculate_(last, nbeads, ljParams_, elParams_, bc_);
        results.push_back(result);
        
        // Collect potential energies and forces.
        std::vector<force_t> forces(nbeads, force_t{});
        energy_t epot{0.0};
        for (auto result : results) {
            const std::vector<force_t>& ff = result.second;
            for (std::size_t index = 0; index != ff.size(); ++index) {
                forces[index] += ff[index];
            }
            epot += result.first;
        }
        
        // Set forces of particles.
        for (std::size_t index = 0; index != forces.size(); ++index) {
            Bead& bead = *all[index];
            bead.force(forces[index]);
        }
        
        return epot;
    }
    
    std::string LJCoulombForces<Bead>::id() const
    {
        return "lj-coulomb-forces";
    }
}
