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
 * File:   model-factory.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 30, 2019, 12:51 PM
 */

#include "simploce/simulation/sim-model-factory.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/pt-langevin-velocity-verlet.hpp"
#include "simploce/simulation/cg-pol-water.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/particle/particle-model-factory.hpp"
#include "simploce/util/mu-units.hpp"
#include <algorithm>
#include <random>
#include <cmath>
#include <memory>

namespace simploce {
    
    const static length_t PT_RMAX = 0.4; // nm.
    
    SimulationModelFactory::SimulationModelFactory(const particle_model_fact_ptr_t& particleModelFactory,
                                                   const spec_catalog_ptr_t& catalog) :
        particleModelFactory_{particleModelFactory}, catalog_{catalog}
    {        
    }
        
    cg_sim_model_ptr_t 
    SimulationModelFactory::polarizableWater(const box_ptr_t& box,
                                             const density_t atDensitySI,
                                             const temperature_t temperature)
    {
        std::clog << "Creating coarse grained simulation model for polarizable water." << std::endl;
        
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        cg_ptr_t cg = 
            particleModelFactory_->polarizableWater(box, atDensitySI, temperature);
        
        // Periodic boundary conditions.
        bc_ptr_t bc = factory::pbc(box);
        std::clog << "Using periodic boundary conditions." << std::endl;

        // Interactor.
        cg_interactor_ptr_t interactor = 
                factory::polarizableWaterInteractor(catalog_, box, bc);
        
        // Displacer.
        std::shared_ptr<CoarseGrainedDisplacer> displacer = 
            std::make_shared<LangevinVelocityVerlet<CoarseGrained>>(interactor);
        std::clog << "Using \"Langevin Velocity Verlet\" algorithm." << std::endl;
        
        // Done.
        return std::make_shared<cg_sim_model_t>(cg, displacer, interactor, box, bc);
    }    
    
    cg_sim_model_ptr_t 
    SimulationModelFactory::formicAcidSolution(const box_ptr_t& box,
                                               const density_t atDensitySI,
                                               const molarity_t molarity,
                                               const temperature_t temperature,
                                               bool protonatable)
    {
        std::clog << "Creating coarse grained simulation model for formic acid (HCOOH)." 
                  << std::endl;

        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        std::clog << "HCOOH Molarity: " << molarity << " M" << std::endl;
        
        // Create formic acid in CG polarizable water.
        cg_ptr_t cg = particleModelFactory_->formicAcidSolution(box, 
                                                                atDensitySI, 
                                                                molarity, 
                                                                temperature,
                                                                protonatable);
        
        // Periodic boundary conditions.
        bc_ptr_t bc = factory::pbc(box);
        std::clog << "Using periodic boundary conditions." << std::endl;

        // Interactor.
        cg_interactor_ptr_t interactor =
            factory::formicAcidSolutionInteractor(catalog_, box, bc);            
        
        // Pair list for protonatables.
        pt_pair_list_gen_ptr_t generator = 
            factory::protonTransferPairListGenerator(PT_RMAX, bc);
        
        // Proton transfer.
        rate_t rate = 1.0/1.5;
        real_t gamma = 1.0 / rate();
        pt_displacer_ptr_t ptDisplacer = factory::constantRate(rate, gamma);
        
        // Displacer.
        std::shared_ptr<CoarseGrainedDisplacer> displacer = 
            std::make_shared<ProtonTransferLangevinVelocityVerlet>(interactor, 
                                                                   generator,
                                                                   ptDisplacer);
        std::clog << "Using \"Langevin Velocity Verlet + Proton Transfer\" algorithm." << std::endl;
        
        // Done.
        return std::make_shared<cg_sim_model_t>(cg, displacer, interactor, box, bc);

    }
    
    cg_sim_model_ptr_t SimulationModelFactory::readCoarseGrainedFrom(std::istream& stream)
    {
        cg_sim_model_ptr_t cg{new SimulationModel<Bead>};
        cg->readFrom(stream, catalog_);
        return cg;
    }
  
}
