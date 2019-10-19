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
 * File:   acid-base-solution.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 8, 2019, 3:47 PM
 */

#include "simploce/simulation/acid-base-solution.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/simulation/lj-coulomb-forces.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/util/map2.hpp"
#include "simploce/util/map.hpp"
#include <memory>
#include <stdexcept>
#include <cmath>
#include <map>

namespace simploce {
    
    // Martini force field. Marrink, J. Phys. B. 111, 7812-7824, 2007.
    // The following refers to level I interactions between polar groups in water.
    static real_t EPS = 5.0;                              // kJ/mol
    static real_t SIGMA = 0.62;                           // nm.
    static real_t C12 = 4.0 * EPS * std::pow(SIGMA, 12);  // kJ nm^12/ mol
    static real_t C6 = 4.0 * EPS * std::pow(SIGMA, 6);    // kJ nm^6 / mol
    
    using lj_params_t = ForceField::lj_params_t;
    using el_params_t = ForceField::el_params_t;
    
    static lj_params_t ljParams_;
    static el_params_t elParams_;
    
    static std::unique_ptr<LJCoulombForces<Bead>> LJ_COULOMB_F{};
    
    void
    setup_(const spec_catalog_ptr_t& catalog,
           const bc_ptr_t& bc,
           const cg_ff_ptr_t& water)
    {
        // Polarizable water.
        spec_ptr_t PCW = catalog->lookup("PCW");
        spec_ptr_t DP = catalog->lookup("DP");
        
        // Acids.
        spec_ptr_t HCOOH = catalog->lookup("HCOOH");
        
        // Bases.
        spec_ptr_t NH4 = catalog->lookup("NH4");
        
        // Water parameters.
        auto parameters = water->parameters();
        ljParams_ = parameters.first;
        elParams_ = parameters.second;
        
        // LJ, 
        //auto zero = std::make_pair(0.0, 0.0);    // No/zero interaction parameters.
        
        // water: correct for protonatable.
        if ( !ljParams_.contains("PCW", "PCW") ) {
            throw std::domain_error(
                "Missing LJ parameters for polarizable/protontable water."
            );
        }
        auto PCW_PCW = ljParams_.at("PCW", "PCW");
        //ljParams_.add(PCW->name(), DP->name(), zero);
        //ljParams_.add(DP->name(), PCW->name(), zero);
        //ljParams_.add(DP->name(), DP->name(), zero);
        
        auto c12 = PCW_PCW.first;
        auto c6 = PCW_PCW.second;
        real_t sigma = std::pow(c12/c6, 1.0/6.0);
        real_t eps = c6 * c6 / (4.0 * c12);
        
        // HCOOH-HCOOH, HCOOH-water
        auto HCOOH_HCOOH = std::make_pair(C12, C6);
        ljParams_.add(HCOOH->name(), PCW->name(), HCOOH_HCOOH);
        auto s = (sigma + SIGMA) / 2.0;
        auto e = std::sqrt(eps + EPS);
        c12 = 4.0 * e * std::pow(s, 12);
        c6 = 4.0 * e * std::pow(s, 6);
        auto HCOOH_PCW = std::make_pair(c12, c6);
        ljParams_.add(PCW->name(), HCOOH->name(), HCOOH_PCW);
        ljParams_.add(HCOOH->name(), PCW->name(), HCOOH_PCW);
        //ljParams_.add(HCOOH->name(), DP->name(), zero);
        //ljParams_.add(DP->name(), HCOOH->name(), zero);
        
        std::clog << "Acids/Bases in polarizable water:" << std::endl;
        std::clog << "Electrostatic interaction parameters:" << std::endl;
        std::clog << elParams_ << std::endl;
        std::clog << "LJ Interaction parameters" << std::endl;
        std::clog << ljParams_ << std::endl;
        
        LJ_COULOMB_F = std::make_unique<LJCoulombForces<Bead>>(ljParams_, elParams_, bc);        
    }
    
    AcidBaseSolution::AcidBaseSolution(const spec_catalog_ptr_t& catalog,
                                       const bc_ptr_t& bc,
                                       const cg_ff_ptr_t& water) :
        catalog_{catalog}, bc_{bc}, water_{water}
    {   
        if ( !catalog_ ) {
            throw std::domain_error(
                "AcidBaseSolution: Missing particle specification catalog."
            );
        }
        if ( !bc_ ) {
            throw std::domain_error(
                "AcidBaseSolution: Missing boundary condition."
            );            
        }
        if ( !water ) {
            throw std::domain_error(
                "AcidBaseSolution: Missing protonatable water."
            );                        
        }
        
        setup_(catalog_, bc_, water_);
    }
    
    energy_t 
    AcidBaseSolution::interact(const std::vector<bead_ptr_t>& all,
                               const std::vector<bead_ptr_t>& free,
                               const std::vector<bead_group_ptr_t>& groups,
                               const std::vector<bead_pair_list_t>& pairLists)
    {
        energy_t epot = water_->bonded(all, free, groups, pairLists);
        epot += LJ_COULOMB_F->interact(all, free, groups, pairLists);
        return epot;        
    }
    
    energy_t 
    AcidBaseSolution::interact(const bead_ptr_t& bead,
                               const std::vector<bead_ptr_t>& all,
                               const std::vector<bead_ptr_t>& free,
                               const std::vector<bead_group_ptr_t>& groups)
    {
        return energy_t{};
    }
    
    energy_t 
    AcidBaseSolution::bonded(const std::vector<bead_ptr_t>& all,
                             const std::vector<bead_ptr_t>& free,
                             const std::vector<bead_group_ptr_t>& groups,
                             const std::vector<bead_pair_list_t>& pairLists)
    {
        return water_->bonded(all, free, groups, pairLists);
    }
    
    std::string 
    AcidBaseSolution::id() const
    {
        return conf::ACID_BASE_SOLUTION;
    }
    
    std::pair<lj_params_t, el_params_t> 
    AcidBaseSolution::parameters() const
    {
        return water_->parameters();
    }
}
