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
#include <utility>
#include <memory>

namespace simploce {
    
    // Interaction parameters, from Riniker et al, 2011.
    static const real_t EPS_R = 78.5;             // Relative permittivity
    static const length_t R_CW_DP = 0.2;          // nm.
    static particle_spec_ptr_t CW{};
    static particle_spec_ptr_t DP{};
    static const real_t FC = 2.0e+06;             // Force constant in kJ/(mol nm^4)
    static const real_t C12_CW_CW = 1.298e-03;    // kJ/(mol nm^12)
    static const real_t C6_CW_CW = 0.088;         // kJ/(mol nm^6)
    
    static std::unique_ptr<LJCoulombForces<Bead>> LJ_COULOMB_F{};
    
    static void setup_(const spec_catalog_ptr_t& catalog,
                       const bc_ptr_t& bc)
    {
        CW = catalog->lookup("CW");
        DP = catalog->lookup("DP");
    
        // Electrostatics.
        LJCoulombForces<Bead>::el_params_t elParams;
        auto eps_r = std::make_pair("eps_r", EPS_R);
        elParams.insert(eps_r);

        // LJ
        LJCoulombForces<Bead>::lj_params_t ljParams;
        auto CW_CW = std::make_pair(C12_CW_CW, C6_CW_CW);    
        ljParams.add(CW->name(), CW->name(), CW_CW);
        auto zero = std::make_pair(0.0, 0.0);    
        ljParams.add(CW->name(), DP->name(), zero);
        ljParams.add(DP->name(), CW->name(), zero);
        ljParams.add(DP->name(), DP->name(), zero);

        std::clog << "LJ Interaction parameters" << std::endl;
        std::clog << ljParams << std::endl;
    
        LJ_COULOMB_F = std::make_unique<LJCoulombForces<Bead>>(ljParams, elParams, bc);
  }

    
    static std::pair<energy_t, std::vector<force_t>> 
    bonded_(std::size_t nparticles, const std::vector<bead_group_ptr_t>& groups)
    {
        using bond_cont_t = typename ParticleGroup<Bead>::bond_cont_t;
        
        static real_t fc2 = 2.0 * FC;
        static real_t halve_fc = 0.5 * FC;
        static force_t fi;
        static force_t fj;
        
        // Potential energy and forces.
        energy_t epot{0.0};
        std::vector<force_t> forces(nparticles, force_t{});
        
        for (auto g : groups) {
            ParticleGroup<Bead>& group = *g;
            const bond_cont_t& bonds = group.bonds();
            const Bond<Bead>& bond = bonds[0];                    // There is 1 bond only.
            const bead_ptr_t pi = bond.getParticleOne();
            std::size_t index_i = pi->id() - 1;                   // Particle id starts at 1.
            const bead_ptr_t pj = bond.getParticleTwo();
            std::size_t index_j = pj->id() - 1;                   // Particle id starts at 1.
      
            position_t ri = pi->position();
            position_t rj = pj->position();
            dist_vect_t rij = ri - rj;            // Distance vector (nm), no boundary conditions.
            length_t R = norm<length_t>(rij);     // Distance (nm)

            // Unit distance vector.
            dist_vect_t uv = rij / R;

            length_t dis = R - R_CW_DP;
            if ( dis() > 0.0 ) {
                real_t dis3 = dis() * dis() * dis();
                real_t dis4 = dis() * dis3; 
                epot += halve_fc * dis4;
                real_t f = -fc2 * dis3;
                for (std::size_t k = 0; k != 3; ++k) {
                    real_t fv = f * uv[k];
                    fi[k] = fv;
                    fj[k] = -fv;
                }
                forces[index_i] += fi;
                forces[index_j] += fj;
            }
        }
        return std::make_pair(epot, forces);
    }
    
    CoarseGrainedPolarizableWater::CoarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog,
                                                                 const bc_ptr_t& bc) :
        CoarseGrainedForceField{}, catalog_{catalog}, bc_{bc}
    {      
            setup_(catalog, bc);
    }
    
    energy_t CoarseGrainedPolarizableWater::interact(const std::vector<bead_ptr_t>& all,
                                                     const std::vector<bead_ptr_t>& free,
                                                     const std::vector<bead_group_ptr_t>& groups,
                                                     const std::vector<bead_pair_list_t>& pairLists)
    {
        std::pair<energy_t, std::vector<force_t>> bonded = bonded_(all.size(), groups);
        energy_t epot = LJ_COULOMB_F->interact(all, free, groups, pairLists);
        epot += bonded.first;
        
        std::vector<force_t>& forces = bonded.second;
        for (std::size_t index = 0; index != all.size(); ++index) {
            Bead& bead = *all[index];
            force_t force = bead.force();
            force += forces[index];
            bead.force(force);
        }
        
        return epot;
    }
    
    std::string CoarseGrainedPolarizableWater::id() const
    {
        return conf::POLARIZABLE_WATER;
    }
    
    length_t CoarseGrainedPolarizableWater::idealDistanceCWDP()
    {
        return R_CW_DP;
    }
}
