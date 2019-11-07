/*
 * The MIT License
 *
 * Copyright 2019 juffer.
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

#include "simploce/simulation/cg-lj-fluid.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/simulation/lj-coulomb-forces.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/map2.hpp"
#include "simploce/util/map.hpp"
#include "simploce/util/mu-units.hpp"

namespace simploce {
    
    using lj_params_t = ForceField::lj_params_t;
    using el_params_t = ForceField::el_params_t;
    
    static const real_t KB = MUUnits<real_t>::KB;
    
    // Values are taken from Marrink et al, J. Phys. Chem. B 2007, 111, 7812-7824
    // used the nonpolar interactions between beads in aliphatic chains
    static const real_t EPS_R = 78.5;               // Relative permittivity.
    static const energy_t EPS = 3.5;                // LJ eps in kJ/mol.
    static const length_t SIGMA = 0.47;             // LJ sigma in nm.
    
    static std::unique_ptr<LJCoulombForces<Bead>> LJ_COULOMB_F{};
    
    static el_params_t elParams_;
    static lj_params_t ljParams_;
    
    static void setup_(const spec_catalog_ptr_t& catalog,
                       const bc_ptr_t& bc,
                       const box_ptr_t& box)
    {
        // Electrostatics. Not used, but must be present.
        auto eps_r = std::make_pair("eps_r", EPS_R);
        elParams_.insert(eps_r);

        // LJ.
        spec_ptr_t spec = catalog->lookup("AP");

        real_t C12 = 4.0 * EPS() * std::pow(SIGMA(), 12.0);
        real_t C6 = 4.0 * EPS() * std::pow(SIGMA(), 6.0);
        auto pair = std::make_pair(C12, C6);
        ljParams_.add(spec->name(), spec->name(), pair);

        std::clog << "LJ Fluid:" << std::endl;
        std::clog << "LJ Interaction parameters" << std::endl;
        std::clog << ljParams_ << std::endl;        

        LJ_COULOMB_F = 
            std::make_unique<LJCoulombForces<Bead>>(ljParams_, elParams_, bc, box);
    }
    
    CoarseGrainedLJFluid::CoarseGrainedLJFluid(const spec_catalog_ptr_t& catalog,
                                               const bc_ptr_t& bc,
                                               const box_ptr_t& box) :
        catalog_{catalog}, bc_{bc}, box_{box}
    {        
            setup_(catalog_, bc_, box_);
    }
    
    std::pair<energy_t, energy_t> 
    CoarseGrainedLJFluid::interact(const std::vector<bead_ptr_t>& all,
                                   const std::vector<bead_ptr_t>& free,
                                   const std::vector<bead_group_ptr_t>& groups,
                                   const PairLists<Bead>& pairLists)
    {
        return LJ_COULOMB_F->interact(all, free, groups, pairLists);
    }
    
    energy_t 
    CoarseGrainedLJFluid::bonded(const std::vector<bead_ptr_t>& all,
                                 const std::vector<bead_ptr_t>& free,
                                 const std::vector<bead_group_ptr_t>& groups,
                                 const PairLists<Bead>& pairLists)
    {
        // There are no bonded interaction.
        return energy_t{0.0};
    }
    
    std::pair<energy_t, energy_t> 
    CoarseGrainedLJFluid::interact(const bead_ptr_t& bead,
                                   const std::vector<bead_ptr_t>& all,
                                   const std::vector<bead_ptr_t>& free,
                                   const std::vector<bead_group_ptr_t>& groups) 
    {
        return LJ_COULOMB_F->interact(bead, all, free, groups);
    }
    
    std::string 
    CoarseGrainedLJFluid::id() const
    {
        return conf::LJ_FLUID;
    }
    
    std::pair<lj_params_t, el_params_t> 
    CoarseGrainedLJFluid::parameters() const
    {
        return std::pair<lj_params_t, el_params_t>{};
    }
}