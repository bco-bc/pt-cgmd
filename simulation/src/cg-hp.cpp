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

#include "simploce/simulation/cg-hp.hpp"
#include "simploce/simulation/sconf.hpp"
#include <tuple>

namespace simploce {
    
    using lj_params_t = ForceField::lj_params_t;
    using el_params_t = ForceField::el_params_t;
    using result_t = std::pair<energy_t, std::vector<force_t>>;
    
    static const real_t FC = 1000.0;    
    static const length_t R_REF = 0.4;  // nm.
    
    static energy_t
    bonded_(const bead_group_ptr_t& group,
            std::vector<force_t>& forces)
    {
        using bond_cont_t = typename ParticleGroup<Bead>::bond_cont_t;
        
        energy_t epot = 0;
        const bond_cont_t& bonds = group->bonds();
        for (auto& bond : bonds) {
        
            // First particle.
            const auto pi = bond.getParticleOne();
            std::size_t index_i = pi->index();
            position_t ri = pi->position();
        
            // Second particle.
            const auto pj = bond.getParticleTwo();
            std::size_t index_j = pj->index();      
            position_t rj = pj->position();

            dist_vect_t rij = ri - rj;
            real_t Rij = norm<real_t>(rij);
        
            real_t R = Rij - R_REF();
        
            epot += 0.5 * FC * R * R;
            dist_vect_t uv = rij/Rij;                          // Unit vector
            real_t dHPdR = FC * R;
            force_t f{};
            for (std::size_t k = 0; k != 3; ++k) {
                f[k] = -dHPdR * uv[k];
            }
        
            // Store energy and forces.
            forces[index_i] += f;
            forces[index_j] -= f;
        }
        
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
    
    HarmonicPotential::HarmonicPotential(const spec_catalog_ptr_t& catalog,
                                         const bc_ptr_t& bc,
                                         const box_ptr_t& box) : 
        catalog_{catalog}, bc_{bc}, box_{box}
    {            
    }
        
    energy_t HarmonicPotential::interact(const std::vector<bead_ptr_t>& all,
                                         const std::vector<bead_ptr_t>& free,
                                         const std::vector<bead_group_ptr_t>& groups,
                                         const PairLists<Bead>& pairLists)
    {
        return this->bonded(all, free, groups, pairLists);
    }
    
    energy_t HarmonicPotential::bonded(const std::vector<bead_ptr_t>& all,
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
    
    energy_t HarmonicPotential::interact(const bead_ptr_t& bead,
                                         const std::vector<bead_ptr_t>& all,
                                         const std::vector<bead_ptr_t>& free,
                                         const std::vector<bead_group_ptr_t>& groups)
    {
        // Forces are not used.
        std::vector<force_t> forces(all.size(), force_t{});
        
        energy_t epot{0.0};
        for (auto g : groups) {            
            if ( g->contains(bead) ) {
                epot += bonded_(g, forces);
            }
        }
        return epot;
    }
    
    std::string HarmonicPotential::id() const
    {
        return conf::HP;
    }
        
    std::pair<lj_params_t, el_params_t> HarmonicPotential::parameters() const
    {
        return std::pair<lj_params_t, el_params_t>{};
    }
    
}