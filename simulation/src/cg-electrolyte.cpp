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
 * File:   cg-electrolyte.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 18, 2019, 1:04 PM
 */

#include "simploce/simulation/cg-electrolyte.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/simulation/lj-coulomb-forces.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/util/map2.hpp"
#include "simploce/util/map.hpp"
#include "simploce/util/mu-units.hpp"

namespace simploce {
    
    using lj_params_t = ForceField::lj_params_t;
    using el_params_t = ForceField::el_params_t;
    
    static const real_t KB = MUUnits<real_t>::KB;

    // Values are taken from table I, Lenart et al, 2007.
    static const real_t EPS_R = 78.5;                     // Relative permittivity
    static const energy_t NA_NA_EPS = 14.3288 * KB;       // LJ eps in kJ/mol.
    static const length_t NA_NA_SIGMA = 2.443 * 0.1;      // LJ sigma in nm.
    static const energy_t NA_CL_EPS = 42.4104 * KB;       // LJ eps in kJ/mol.
    static const length_t NA_CL_SIGMA = 2.796 * 0.1;      // LJ sigma in nm.
    static const energy_t CL_CL_EPS = 117.7604 * KB;      // LJ eps in kJ/mol.
    static const length_t CL_CL_SIGMA = 3.487 * 0.1;      // LJ sigma in nm.

    static std::unique_ptr<LJCoulombForces<Bead>> LJ_COULOMB_F{};
    
    static el_params_t elParams_;
    static lj_params_t ljParams_;
    
    static void setup_(const spec_catalog_ptr_t& catalog,
                       const bc_ptr_t& bc,
                       const box_ptr_t& box)
    {
        // Electrostatics.
        auto eps_r = std::make_pair("eps_r", EPS_R);
        elParams_.insert(eps_r);

        // LJ.
        spec_ptr_t Na = catalog->lookup("Na+");
        spec_ptr_t Cl = catalog->lookup("Cl-");

        real_t C12_Na_Na = 4.0 * NA_NA_EPS() * std::pow(NA_NA_SIGMA(), 12.0);
        real_t C6_Na_Na = 4.0 * NA_NA_EPS() * std::pow(NA_NA_SIGMA(), 6.0);
        auto Na_Na = std::make_pair(C12_Na_Na, C6_Na_Na);
        ljParams_.add(Na->name(), Na->name(), Na_Na);

        real_t C12_Cl_Cl = 4.0 * CL_CL_EPS() * std::pow(CL_CL_SIGMA(), 12.0);
        real_t C6_Cl_Cl = 4.0 * CL_CL_EPS() * std::pow(CL_CL_SIGMA(), 6.0);
        auto Cl_Cl = std::make_pair(C12_Cl_Cl, C6_Cl_Cl);
        ljParams_.add(Cl->name(), Cl->name(), Cl_Cl);
    
        real_t C12_Na_Cl = 4.0 * NA_CL_EPS() * std::pow(NA_CL_SIGMA(), 12.0);
        real_t C6_Na_Cl = 4.0 * NA_CL_EPS() * std::pow(NA_CL_SIGMA(), 6.0);
        auto Na_Cl = std::make_pair(C12_Na_Cl, C6_Na_Cl);
        ljParams_.add(Na->name(), Cl->name(), Na_Cl);
        ljParams_.add(Cl->name(), Na->name(), Na_Cl);
        
        std::clog << "Electrolyte:" << std::endl;
        std::clog << "Electrostatic interaction parameters:" << std::endl;
        std::clog << elParams_ << std::endl;
        std::clog << "LJ Interaction parameters" << std::endl;
        std::clog << ljParams_ << std::endl;        

        LJ_COULOMB_F = 
            std::make_unique<LJCoulombForces<Bead>>(ljParams_, elParams_, bc, box);
    }    
    
    
    CoarseGrainedElectrolyte::CoarseGrainedElectrolyte(const spec_catalog_ptr_t& catalog,
                                                       const bc_ptr_t& bc,
                                                       const box_ptr_t& box) :
        catalog_{catalog}, bc_{bc}, box_{box}
    {        
            setup_(catalog_, bc_, box_);
    }
    
    energy_t 
    CoarseGrainedElectrolyte::interact(const std::vector<bead_ptr_t>& all,
                                       const std::vector<bead_ptr_t>& free,
                                       const std::vector<bead_group_ptr_t>& groups,
                                       const PairLists<Bead>& pairLists)
    {
        return LJ_COULOMB_F->interact(all, free, groups, pairLists);
    }
    
    energy_t 
    CoarseGrainedElectrolyte::bonded(const std::vector<bead_ptr_t>& all,
                                     const std::vector<bead_ptr_t>& free,
                                     const std::vector<bead_group_ptr_t>& groups,
                                     const PairLists<Bead>& pairLists)
    {
        // There are no bounded interaction.
        return energy_t{0.0};
    }
    
    energy_t 
    CoarseGrainedElectrolyte::interact(const bead_ptr_t& bead,
                                       const std::vector<bead_ptr_t>& all,
                                       const std::vector<bead_ptr_t>& free,
                                       const std::vector<bead_group_ptr_t>& groups) 
    {
        return LJ_COULOMB_F->interact(bead, all, free, groups);
    }
    
    std::string 
    CoarseGrainedElectrolyte::id() const
    {
        return conf::ELECTROLYTE;
    }
    
    std::pair<lj_params_t, el_params_t> 
    CoarseGrainedElectrolyte::parameters() const
    {
        return std::pair<lj_params_t, el_params_t>{};
    }
}