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
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/mu-units.hpp"
#include <future>

namespace simploce {
    
    using lj_params_t = LJCoulombForces<Bead>::lj_params_t;
    using el_params_t = LJCoulombForces<Bead>::el_params_t;
    using result_t = std::pair<energy_t, std::vector<force_t>>;
    using bead_pair_list_t = ParticlePairListGenerator<Bead>::p_pair_list_t;
    
    /**
     * Returns potential energy and force on particle i.
     */
    static std::pair<energy_t, force_t> ljCoulombForce_(const position_t& ri,
                                                        const charge_t& qi,
                                                        const position_t& rj,
                                                        const charge_t& qj,
                                                        real_t C12,
                                                        real_t C6,
                                                        real_t eps_r,
                                                        const bc_ptr_t& bc)
    {
        const real_t four_pi_e0 = MUUnits<real_t>::FOUR_PI_E0;
        
        // Apply boundary condition.
        dist_vect_t rij = bc->apply(ri, rj);
    
        // Potential energy.
        real_t Rij = norm<real_t>(rij);
        real_t Rij2 = Rij * Rij;
        real_t Rij6 = Rij2 * Rij2 * Rij2;
        real_t Rij12 = Rij6 * Rij6;
        real_t t1 = C12 / Rij12;
        real_t t2 = C6 / Rij6;
        real_t LJ = t1 - t2;                          // kj/mol
        real_t t3 = qi * qj / (four_pi_e0 * eps_r);
        real_t elec = t3 / Rij;                       // kJ/mol
        energy_t total = LJ + elec;
    
        // Forces.
        dist_vect_t uv = rij/Rij;
        real_t dLJdR = -6.0 * ( 2.0 * t1 - t2 ) / Rij;  // kJ/(mol nm)
        real_t dElecdR = -t3 / Rij2;
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            real_t fLJ = -dLJdR * uv[k];                 // kJ/(mol nm)
            real_t fElec = -dElecdR * uv[k];             // kJ/(mol nm)
            f[k] = fLJ + fElec;
        }

        return std::make_pair(total, f);
  }    
    
    static result_t forces_(const bead_pair_list_t& single,
                            std::size_t nbeads,
                            const lj_params_t& ljParams,
                            const el_params_t& elParams,
                            const bc_ptr_t& bc)
    {
        std::vector<force_t> forces(nbeads, force_t{});
        energy_t epot{0.0};
    
        // Electrostatic parameters.
        real_t eps_r = elParams.at("eps_r");
    
        for (auto pp : single) {
            
            // First particle
            bead_ptr_t pi = pp.first;
            position_t ri = pi->position();
            std::string name_i = pi->spec()->name();
            charge_t qi = pi->charge();
            std::size_t index_i = pi->id() - 1; // particle id starts at 1!
      
            // Second particle.
            bead_ptr_t pj = pp.second;
            position_t rj = pj->position();
            std::string name_j = pj->spec()->name();
            charge_t qj = pj->charge();
            std::size_t index_j = pj->id() - 1; // particle id starts at 1!

            // LJ parameters.
            auto ljParam = ljParams.at(name_i, name_j);
            energy_t C12 = ljParam.first;
            length_t C6 = ljParam.second;

            std::pair<energy_t, force_t> ef = 
                ljCoulombForce_(ri, qi, rj, qj, C12(), C6(), eps_r, bc);

            // Store.
            epot += ef.first;
            forces[index_i] += ef.second;
            forces[index_j] -= ef.second;
        }
    
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
            std::size_t ntasks = pairLists.size() - 1;  // One other for the current thread.
            futures.clear();
            results.clear();
            for (std::size_t k = 0; k != ntasks; ++k) {
                const bead_pair_list_t& single = pairLists[k];
                futures.push_back(
                    std::async(
                        std::launch::async, 
                        forces_,
                        std::ref(single),
                        nbeads,
                        std::ref(ljParams_),
                        std::ref(elParams_),
                        std::ref(bc_)
                    )
                );
            }
            // Wait for these tasks to complete.
            results = util::waitForAll<result_t>(futures);
        }
    
        // Current thread.
        const bead_pair_list_t& last = pairLists[nlists - 1];
        result_t result = forces_(last, nbeads, ljParams_, elParams_, bc_);
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
