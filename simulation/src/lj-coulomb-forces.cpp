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
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/simulation/sim-util.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/mu-units.hpp"
#include <future>
#include <utility>
#include <cassert>

namespace simploce {
    
    using lj_params_t = ForceField::lj_params_t;
    using el_params_t = ForceField::el_params_t;
    using result_t = std::pair<energy_t, std::vector<force_t>>;
    using bead_pair_list_t = PairLists<Bead>::pp_list_cont_t;
    using bead_group_pair_list_t = PairLists<Bead>::pg_list_cont_t;
    using group_pair_list_t = PairLists<Bead>::gg_list_cont_t;
    
    /**
     * Returns interaction potential energy and force on particle i.
     */
    static std::pair<energy_t, force_t> ljCoulombForce_(const position_t& ri,
                                                        const charge_t& qi,
                                                        const position_t& rj,
                                                        const charge_t& qj,
                                                        real_t C12,
                                                        real_t C6,
                                                        real_t eps_r,
                                                        const bc_ptr_t& bc,
                                                        bool applyBC,
                                                        const dist_vect_t& R)
    {
        static const real_t four_pi_e0 = MUUnits<real_t>::FOUR_PI_E0;
        

        // Apply boundary condition.
        dist_vect_t rij = applyBC ? bc->apply(ri, rj) : ri - rj - R;
        real_t Rij = norm<real_t>(rij);

        // Potential energy
        
        // LJ
        real_t Rij2 = Rij * Rij;
        real_t Rij6 = Rij2 * Rij2 * Rij2;
        real_t Rij12 = Rij6 * Rij6;
        real_t t1 = C12 / Rij12;
        real_t t2 = C6 / Rij6;
        real_t LJ = t1 - t2;                          // kj/mol
        
        // Coulomb
        real_t t3 = qi * qj / (four_pi_e0 * eps_r);
        real_t elec = t3 / Rij;  // kJ/mol
        
        // Total potential energy.
        energy_t epot = LJ + elec;
    
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

        return std::make_pair(epot, f);
    }
    
    // Returns forces on beads and energy for bead pairs.
    static result_t ppForces_(const bead_pair_list_t ppPairList,
                              std::size_t nbeads,
                              const lj_params_t& ljParams,
                              const el_params_t& elParams,
                              const bc_ptr_t& bc)
    {
        std::vector<force_t> forces(nbeads, force_t{});
        energy_t epot{0.0};
        static const dist_vect_t R{};
        
        // Electrostatic parameters.
        static const real_t eps_r = elParams.at("eps_r");
        
    
        for (auto pp : ppPairList) {
            
            // First particle
            bead_ptr_t pi = pp.first;
            position_t ri = pi->position();
            std::string name_i = pi->spec()->name();
            charge_t qi = pi->charge();
            std::size_t index_i = pi->index();
      
            // Second particle.
            bead_ptr_t pj = pp.second;
            position_t rj = pj->position();
            std::string name_j = pj->spec()->name();
            charge_t qj = pj->charge();
            std::size_t index_j = pj->index();
            
            // Calculate interaction.
            auto ljParam = ljParams.at(name_i, name_j);
            auto C12 = ljParam.first;
            auto C6 = ljParam.second;
            auto ef = ljCoulombForce_(ri, qi, rj, qj, C12, C6, eps_r, bc, true, R);

            // Store energy and forces.
            epot += ef.first;
            forces[index_i] += ef.second;
            forces[index_j] -= ef.second;
        }
    
        return std::make_pair(epot, forces);
    }
    
    // Returns forces on beads and energy for for bead/bead group pairs.
    static result_t pgForces_(const bead_group_pair_list_t& pgPairList,
                              std::size_t nbeads,
                              const lj_params_t& ljParams,
                              const el_params_t& elParams,
                              const bc_ptr_t& bc)
    {
        std::vector<force_t> forces(nbeads, force_t{});
        energy_t epot{0.0};
        static const dist_vect_t R{};
                
        // Electrostatic parameters.
        static const real_t eps_r = elParams.at("eps_r");
            
        for (auto pg : pgPairList) {
            
            // First particle
            bead_ptr_t pi = pg.first;
            position_t ri = pi->position();
            std::string name_i = pi->spec()->name();
            charge_t qi = pi->charge();
            std::size_t index_i = pi->index();
            
            bead_group_ptr_t g = pg.second;
            for (auto pj : g->particles()) {
            
                // Second particle.
                position_t rj = pj->position();
                std::string name_j = pj->spec()->name();
                charge_t qj = pj->charge();
                std::size_t index_j = pj->index();       
                
                // Calculate interaction.
                auto ljParam = ljParams.at(name_i, name_j);
                auto C12 = ljParam.first;
                auto C6 = ljParam.second;
                auto ef = ljCoulombForce_(ri, qi, rj, qj, C12, C6, eps_r, bc, true, R);

                // Store energy and forces.
                epot += ef.first;
                forces[index_i] += ef.second;
                forces[index_j] -= ef.second;
            }            
        }
        
        return std::make_pair(epot, forces);
    }
    
    // Returns forces on beads and energy for bead group/bead group pairs.
    static result_t ggForces_(const group_pair_list_t& ggPairList,
                              std::size_t nbeads,
                              const lj_params_t& ljParams,
                              const el_params_t& elParams,
                              const bc_ptr_t& bc)
    {
        std::vector<force_t> forces(nbeads, force_t{});
        energy_t epot{0.0};
        
        // Electrostatic parameters.
        static real_t eps_r = elParams.at("eps_r");

        for (auto gg : ggPairList) {
            
            // Groups.
            bead_group_ptr_t gi = gg.first;
            position_t rg_i = gi->position();
            bead_group_ptr_t gj = gg.second;
            position_t rg_j = gj->position();
            dist_vect_t Rbc = bc->apply(rg_i, rg_j);
            dist_vect_t R = Rbc - (rg_i - rg_j);
            
            for (auto pi : gi->particles()) {
                
                // First particle
                position_t ri = pi->position();
                std::string name_i = pi->spec()->name();
                charge_t qi = pi->charge();
                std::size_t index_i = pi->index();
                
                for (auto pj : gj->particles()) {
                
                    // Second particle.
                    position_t rj = pj->position();
                    std::string name_j = pj->spec()->name();
                    charge_t qj = pj->charge();
                    std::size_t index_j = pj->index();       
                
                    // Calculate interaction.
                    auto ljParam = ljParams.at(name_i, name_j);
                    auto C12 = ljParam.first;
                    auto C6 = ljParam.second;
                    auto ef = ljCoulombForce_(ri, qi, rj, qj, C12, C6, eps_r, bc, true, R);

                    // Store energy and forces.
                    epot += ef.first;
                    forces[index_i] += ef.second;
                    forces[index_j] -= ef.second;                                        
                }
            }
        }
        
        return std::make_pair(epot, forces);
    }
    
    // Interaction energy only, forces are ignored.
    static energy_t energy_(const bead_ptr_t& bead,
                            const std::vector<bead_ptr_t>& free,
                            const lj_params_t& ljParams,
                            const el_params_t& elParams,
                            const bc_ptr_t& bc,
                            const box_ptr_t& box)
    {
        // Electrostatic parameters.
        static const real_t eps_r = elParams.at("eps_r");

        length_t rc = util::cutoffDistance(box);
        real_t rc2 = rc * rc;
        energy_t epot{0.0};
        static dist_vect_t R{};
        
        position_t ri = bead->position();
        std::string name_i = bead->spec()->name();
        charge_t qi = bead->charge();
        std::size_t index_i = bead->index();

        for (auto f : free) {
            std::size_t index_j = f->index();
            if ( index_j != index_i ) {
                position_t rj = f->position();
                auto rij = bc->apply(ri, rj);
                auto Rij2 = norm2<real_t>(rij);
                if ( Rij2 < rc2 ) {
                    std::string name_j = f->spec()->name();
                    charge_t qj = f->charge();
                                  
                    // Calculate interaction.      
                    auto ljParam = ljParams.at(name_i, name_j);
                    auto C12 = ljParam.first;
                    auto C6 = ljParam.second;
                    auto ef = ljCoulombForce_(ri, qi, rj, qj, C12, C6, eps_r, bc, true, R);
                    
                    // Store interaction energy.
                    epot += ef.first;
                }                
            }
        }
        
        return epot;
    }
    
    energy_t energy_(const bead_ptr_t& bead,
                     const std::vector<bead_group_ptr_t>& groups,
                     const lj_params_t& ljParams,
                     const el_params_t& elParams,
                     const bc_ptr_t& bc,
                     const box_ptr_t& box)
    {
        // Electrostatic parameters.
        static const real_t eps_r = elParams.at("eps_r");

        length_t rc = util::cutoffDistance(box);
        real_t rc2 = rc * rc;
        energy_t epot{0.0};
        static dist_vect_t R{};
        
        position_t ri = bead->position();
        std::string name_i = bead->spec()->name();
        charge_t qi = bead->charge();
        
        for (auto g : groups) {
            if ( !g->contains(bead) ) {
                position_t rj = g->position();
                auto rij = bc->apply(ri, rj);
                auto Rij2 = norm2<real_t>(rij);
                if ( Rij2 < rc2 ) {
                    for (auto p : g->particles()) {
                        rj = p->position();
                        std::string name_j = p->spec()->name();
                        charge_t qj = p->charge();
                        
                        // Calculate interaction.
                        auto ljParam = ljParams.at(name_i, name_j);
                        auto C12 = ljParam.first;
                        auto C6 = ljParam.second;
                        auto ef = ljCoulombForce_(ri, qi, rj, qj, C12, C6, eps_r, bc, true, R);
                        
                        // Store interaction.
                        epot += ef.first;                                                
                    }
                }
            }
        }
        
        return epot;        
    }
    
    LJCoulombForces<Bead>::LJCoulombForces(const lj_params_t& ljParams, 
                                           const el_params_t& elParams, 
                                           const bc_ptr_t& bc,
                                           const box_ptr_t& box) :
        ljParams_{ljParams}, elParams_{elParams}, bc_{bc}, box_{box}
    {        
    }
        
    energy_t LJCoulombForces<Bead>::interact(const std::vector<bead_ptr_t>& all,
                                             const std::vector<bead_ptr_t>& free,
                                             const std::vector<bead_group_ptr_t>& groups,
                                             const PairLists<Bead>& pairLists)
    {         
        static std::vector<group_pair_list_t> subPairLists =
            util::makeSubLists(pairLists.groupPairList());
        static bool firstTime = true;
        
        // Holds all force calculation results.
        std::vector<result_t> results{};
        
        auto nbeads = all.size();
        
        // Concurrent calculation only for large number of particles.        
        if ( nbeads > conf::MIN_NUMBER_OF_PARTICLES ) {
            
            // Concurrently.
            
            std::vector<std::future<result_t> > futures{};
                        
            // Handle particle group/particle group interaction concurrently,
            // where one task is executed by the current thread.
            if ( pairLists.isModified() && !firstTime) {
                subPairLists = util::makeSubLists(pairLists.groupPairList());
                firstTime = false;
            }
            std::size_t ntasks = subPairLists.size() - 1;
            if ( ntasks > 1 ) {
                // Setup concurrent force calculations.
                for (std::size_t k = 0; k != ntasks; ++k) {
                    const auto& single = subPairLists[k];
                    futures.push_back(
                        std::async(
                            std::launch::async, 
                            ggForces_,
                            std::ref(single),
                            nbeads,
                            std::ref(ljParams_),
                            std::ref(elParams_),
                            std::ref(bc_)
                        )
                    );
                }
            
                // Wait for tasks to complete.
                results = util::waitForAll<result_t>(futures);
            }
            
            // The rest is handled sequentially by the current thread
            
            // One remaining particle group/particle group interaction.
            const auto& single = subPairLists[ntasks - 1];
            if ( !single.empty() ) {
                auto ggResult = 
                    ggForces_(single, nbeads, ljParams_, elParams_, bc_);
                results.push_back(ggResult);
            }
            
            // Interaction between free particles.
            auto& ppPairList = pairLists.particlePairList();
            if ( !ppPairList.empty() ) {
                auto ppResult = 
                    ppForces_(ppPairList, nbeads, ljParams_, elParams_, bc_);
                results.push_back(ppResult);
            }
            
            // Interaction between free particles and particle groups.
            auto& pgPairList = pairLists.particleGroupPairList();
            if ( !pgPairList.empty() ) {
                auto pgResult =
                    pgForces_(pgPairList, nbeads, ljParams_, elParams_, bc_);
                results.push_back(pgResult);
            }
        
        } else {
                    
            // Sequentially            
            auto ppResult = 
                ppForces_(pairLists.particlePairList(), 
                          nbeads, 
                          ljParams_, 
                          elParams_, 
                          bc_);
            results.push_back(ppResult);
            auto pgResult =
                pgForces_(pairLists.particleGroupPairList(),
                          nbeads,
                          ljParams_, 
                          elParams_, 
                          bc_);
            results.push_back(pgResult);
            auto ggResult =
                ggForces_(pairLists.groupPairList(),
                          nbeads,
                          ljParams_, 
                          elParams_, 
                          bc_);
            results.push_back(ggResult);
        }
        
        // Collect potential energies and forces.
        std::vector<force_t> forces(nbeads, force_t{});
        energy_t epot{0.0};
        for (auto result : results) {
            const std::vector<force_t>& ff = result.second;
            for (std::size_t index = 0; index != nbeads; ++index) {
                forces[index] += ff[index];
            }
            epot += result.first;
        }

        // Update forces of particles.
        for (std::size_t index = 0; index != nbeads; ++index) {
            auto& bead = *all[index];
            force_t f = forces[index] + bead.force();
            bead.force(f);
        }
        
        // Done.
        return epot;
    }
    
    energy_t 
    LJCoulombForces<Bead>::interact(const bead_ptr_t& bead,
                                    const std::vector<bead_ptr_t>& all,
                                    const std::vector<bead_ptr_t>& free,
                                    const std::vector<bead_group_ptr_t>& groups)
    {
        energy_t epot = energy_(bead, free, ljParams_, elParams_, bc_, box_);
        epot += energy_(bead, groups, ljParams_, elParams_, bc_, box_);
        return epot;
    }
    
    energy_t 
    LJCoulombForces<Bead>::bonded(const std::vector<bead_ptr_t>& all,
                                  const std::vector<bead_ptr_t>& free,
                                  const std::vector<bead_group_ptr_t>& groups,
                                  const PairLists<Bead>& pairLists)
    {
        return 0.0;
    }
    
    std::string LJCoulombForces<Bead>::id() const
    {
        return "lj-coulomb-forces";
    }
    
    std::pair<lj_params_t, el_params_t> 
    LJCoulombForces<Bead>::parameters() const
    {
        return std::make_pair(lj_params_t{}, el_params_t{});
    }
}
