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
 * File:   pol-water-force-field.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 30, 2019, 2:17 PM
 */

#include "simploce/simulation/cg-pol-water.hpp"
#include "simploce/simulation/lj-coulomb-forces.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/util/map2.hpp"
#include "simploce/util/map.hpp"
#include <utility>
#include <memory>
#include <map>
#include <iostream>

namespace simploce {
    
    using lj_params_t = ForceField::lj_params_t;
    using el_params_t = ForceField::el_params_t;
    
    // Interaction parameters, adapted from Riniker et al, 2011.
    static const real_t EPS_R = 2.5;              // Relative permittivity
    static const length_t R_CW_DP = 0.2;          // nm.
    static spec_ptr_t CW{};
    static spec_ptr_t DP{};
    static const real_t FC = 2.0e+06;             // Force constant in kJ/(mol nm^4)
    static const real_t C12_CW_CW = 1.298e-03;    // kJ nm^12/mol
    static const real_t C6_CW_CW = 0.088;         // kJ nm^6 /mol
    
    static std::unique_ptr<LJCoulombForces<Bead>> LJ_COULOMB_F{};
    
    static el_params_t elParams_;
    static lj_params_t ljParams_;
    
    static void setup_(const spec_catalog_ptr_t& catalog,
                       const bc_ptr_t& bc,
                       const box_ptr_t& box,
                       bool protonatable)
    {
        CW = protonatable ? catalog->lookup("PCW") : catalog->lookup("CW");
        DP = catalog->lookup("DP");
    
        // Electrostatics.
        auto eps_r = std::make_pair("eps_r", EPS_R);
        elParams_.insert(eps_r);

        // LJ
        auto CW_CW = std::make_pair(C12_CW_CW, C6_CW_CW);    
        ljParams_.add(CW->name(), CW->name(), CW_CW);
        auto zero = std::make_pair(0.0, 0.0);    
        ljParams_.add(CW->name(), DP->name(), zero);
        ljParams_.add(DP->name(), CW->name(), zero);
        ljParams_.add(DP->name(), DP->name(), zero);

        std::clog << "Polarizable water:" << std::endl;
        std::clog << "Electrostatic interaction parameters:" << std::endl;
        std::clog << elParams_ << std::endl;
        std::clog << "LJ Interaction parameters" << std::endl;
        std::clog << ljParams_ << std::endl;
    
        LJ_COULOMB_F = 
            std::make_unique<LJCoulombForces<Bead>>(ljParams_, elParams_, bc, box);
    }
    
    static energy_t
    bonded_(const bead_group_ptr_t& group,
            std::vector<force_t>& forces)
    {
        using bond_cont_t = typename ParticleGroup<Bead>::bond_cont_t;

        static real_t fc2 = 2.0 * FC;
        static real_t halve_fc = 0.5 * FC;
        static force_t fi;
        static force_t fj;
        
        energy_t epot{0.0};
        
        const bond_cont_t& bonds = group->bonds();
        const Bond<Bead>& bond = bonds[0];          // There is 1 bond only.
        const auto pi = bond.getParticleOne();
        std::size_t index_i = pi->index();
        const auto pj = bond.getParticleTwo();
        std::size_t index_j = pj->index();
      
        position_t ri = pi->position();
        position_t rj = pj->position();
        dist_vect_t rij = ri - rj;                 // Distance vector (nm), 
                                                   // no boundary conditions.
        length_t Rij = norm<length_t>(rij);        // Distance (nm)

        length_t dis = Rij - R_CW_DP;
        if ( dis() > 0.0 ) {
            real_t dis3 = dis() * dis() * dis();
            real_t dis4 = dis() * dis3; 
            epot += halve_fc * dis4;
            real_t f = -fc2 * dis3;        
            dist_vect_t uv = rij / Rij;  // Unit distance vector.
            for (std::size_t k = 0; k != 3; ++k) {
                real_t fv = f * uv[k];
                fi[k] = fv;
                fj[k] = -fv;
            }
            forces[index_i] += fi;
            forces[index_j] += fj;
        }        
        
#ifdef _DEBUG        
        if ( epot() > conf::LARGE ) {
            std::clog << "ERROR: High potential energy in water group. " 
                      << "R - R_CW_DP: " << dis << ", "
                      << "Particle indices: " << index_i << ", " << index_j << ", "
                      << "Particle identifiers: " << pi->id() << ", " << pj->id() << ", "
                      << "Energy: " << epot
                      << std::endl;
            throw std::domain_error(
                "High energy polarizable water group displays high energy. See details above."
            );
        }
#endif
        
        return epot;
    }

    
    static std::pair<energy_t, std::vector<force_t>> 
    bonded_(std::size_t nparticles, const std::vector<bead_group_ptr_t>& groups)
    {        
        // Potential energy and forces.
        energy_t epot{0.0};
        std::vector<force_t> forces(nparticles, force_t{});
                
        for (auto g : groups) {            
            epot += bonded_(g, forces);
        }        
        
        return std::make_pair(epot, forces);
    }
    
    CoarseGrainedPolarizableWater::CoarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog,
                                                                 const bc_ptr_t& bc,
                                                                 const box_ptr_t& box,
                                                                 bool protonatable) :
        CoarseGrainedForceField{}, catalog_{catalog}, bc_{bc}, box_{box}
    {      
            setup_(catalog, bc, box, protonatable);
    }
            
    std::pair<energy_t, energy_t> 
    CoarseGrainedPolarizableWater::interact(const std::vector<bead_ptr_t>& all,
                                            const std::vector<bead_ptr_t>& free,
                                            const std::vector<bead_group_ptr_t>& groups,
                                            const PairLists<Bead>& pairLists)
    {
        auto b = bonded_(all.size(), groups);        
        auto bepot = b.first;

        auto nb = LJ_COULOMB_F->interact(all, free, groups, pairLists);
        auto nbepot = nb.second;
        
        auto& forces = b.second;
        for (auto bead : all) {
            auto index = bead->index();
            auto force = bead->force();
            force += forces[index];
            bead->force(force);
        }
        
        return std::make_pair(bepot, nbepot);
    }
    
    energy_t 
    CoarseGrainedPolarizableWater::bonded(const std::vector<bead_ptr_t>& all,
                                          const std::vector<bead_ptr_t>& free,
                                          const std::vector<bead_group_ptr_t>& groups,
                                          const PairLists<Bead>& pairLists)
    {
        auto bonded = bonded_(all.size(), groups);
        auto& forces = bonded.second;
        for (auto bead : all) {
            auto index = bead->index();
            auto force = bead->force();
            force += forces[index];
            bead->force(force);
        }
        
        return bonded.first;
    }
    
    std::pair<energy_t, energy_t>
    CoarseGrainedPolarizableWater::interact(const bead_ptr_t& bead,
                                            const std::vector<bead_ptr_t>& all,
                                            const std::vector<bead_ptr_t>& free,
                                            const std::vector<bead_group_ptr_t>& groups)
    {
        // Forces are not used.
        std::vector<force_t> forces(all.size(), force_t{});
        
        auto nb = LJ_COULOMB_F->interact(bead, all, free, groups);
        energy_t bepot{0.0};
        for (auto g : groups) {            
            if ( g->contains(bead) ) {
                bepot += bonded_(g, forces);
            }
        }
        return std::make_pair(bepot, nb.second);
    }
    
    std::string 
    CoarseGrainedPolarizableWater::id() const
    {
        return conf::POLARIZABLE_WATER;
    }
    
    std::pair<lj_params_t, el_params_t> 
    CoarseGrainedPolarizableWater::parameters() const
    {        
        return std::make_pair(ljParams_, elParams_);
    }
    
    length_t
    CoarseGrainedPolarizableWater::idealDistanceCWDP()
    {
        return R_CW_DP;
    }
}
