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
#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/cg-pol-water.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/particle/particle-model-factory.hpp"
#include "simploce/util/mu-units.hpp"
#include <algorithm>
#include <random>
#include <cmath>
#include <memory>

namespace simploce {
    
    SimulationModelFactory::SimulationModelFactory(const particle_model_fact_ptr_t& particleModelFactory,
                                                   const spec_catalog_ptr_t& catalog) :
        particleModelFactory_{particleModelFactory}, catalog_{catalog}
    {        
    }
        
    cg_sim_model_ptr_t 
    SimulationModelFactory::polarizableWater(const box_ptr_t& box,
                                             const density_t atDensitySI,
                                             const temperature_t temperature,
                                             std::size_t nlimit)
    {
        std::clog << "Creating coarse grained simulation model for polarizable water." << std::endl;
        
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        cg_ptr_t cg = 
            particleModelFactory_->polarizableWater(box, 
                                                    atDensitySI, 
                                                    temperature,
                                                    false,
                                                    nlimit);
        
        // Periodic boundary conditions.
        bc_ptr_t bc = factory::pbc(box);
        std::clog << "Using periodic boundary conditions." << std::endl;

        // Interactor.
        cg_interactor_ptr_t interactor = 
                factory::polarizableWaterInteractor(catalog_, box, bc);
        
        // Displacer.
        std::shared_ptr<CoarseGrainedDisplacer> displacer = 
            //std::make_shared<LeapFrog<CoarseGrained>>(interactor);
        //std::clog << "Using \"Leapfrog\" algorithm." << std::endl;
            std::make_shared<LangevinVelocityVerlet<CoarseGrained>>(interactor);
        std::clog << "Using \"Langevin Velocity Verlet\" algorithm." << std::endl;
        
        // Done.
        return std::make_shared<cg_sim_model_t>(cg, displacer, interactor, box, bc);
    }    
    
    cg_sim_model_ptr_t 
    SimulationModelFactory::formicAcidSolution(const box_ptr_t& box,
                                               bool protonatable,
                                               const molarity_t molarity,
                                               const density_t atDensitySI,
                                               const temperature_t temperature)
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
        
        // Pair list generator for protonatables.
        pt_pair_list_gen_ptr_t generator = 
            factory::protonTransferPairListGenerator(bc);
        
        // Proton transfer.
        pt_displacer_ptr_t ptDisplacer = factory::protonTransferDisplacer();
        
        // Displacer.
        std::shared_ptr<CoarseGrainedDisplacer> displacer = 
            std::make_shared<ProtonTransferLangevinVelocityVerlet>(interactor, 
                                                                   generator,
                                                                   ptDisplacer);
        std::clog << "Using \"Langevin Velocity Verlet + Proton Transfer\" algorithm." << std::endl;
        
        // Done.
        return std::make_shared<cg_sim_model_t>(cg, displacer, interactor, box, bc);

    }
    
    cg_sim_model_ptr_t 
    SimulationModelFactory::electrolyte(const box_ptr_t& box,
                                        molarity_t molarity,
                                        temperature_t temperature)
    {
        std::clog << "Creating coarse grained simulation model for an"
                     " electrolyte solution." << std::endl;
        
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        cg_ptr_t cg = 
            particleModelFactory_->electrolyte(box, molarity, temperature);
        
        // Periodic boundary conditions.
        bc_ptr_t bc = factory::pbc(box);
        std::clog << "Using periodic boundary conditions." << std::endl;

        // Interactor.
        cg_interactor_ptr_t interactor = 
            factory::electrolyteInteractor(catalog_, box, bc);
        
        // Displacer.
        std::shared_ptr<CoarseGrainedDisplacer> displacer = 
            std::make_shared<LeapFrog<CoarseGrained>>(interactor);
        std::clog << "Using \"Leapfrog\" algorithm." << std::endl;
            //std::make_shared<LangevinVelocityVerlet<CoarseGrained>>(interactor);
        //std::clog << "Using \"Langevin Velocity Verlet\" algorithm." << std::endl;
        
        // Done.
        return std::make_shared<cg_sim_model_t>(cg, displacer, interactor, box, bc);                
    }
    
    cg_sim_model_ptr_t 
    SimulationModelFactory::ljFluid(const box_ptr_t& box,
                                    const density_t densitySI,
                                    const temperature_t temperature)
    {
        std::clog << "Creating coarse grained simulation model for a"
                  << " LJ fluid." << std::endl;
        
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        cg_ptr_t cg = 
            particleModelFactory_->ljFluid(box, densitySI, temperature);
            
        // Periodic boundary conditions.
        bc_ptr_t bc = factory::pbc(box);
        std::clog << "Using periodic boundary conditions." << std::endl;

        // Interactor.
        cg_interactor_ptr_t interactor = 
            factory::ljFluidInteractor(catalog_, box, bc);
        
        // Displacer.
        std::shared_ptr<CoarseGrainedDisplacer> displacer = 
            std::make_shared<LeapFrog<CoarseGrained>>(interactor);
                
        // Done.
        return std::make_shared<cg_sim_model_t>(cg, displacer, interactor, box, bc);
    }
    
    cg_sim_model_ptr_t 
    SimulationModelFactory::readCoarseGrainedFrom(std::istream& stream)
    {
        cg_sim_model_ptr_t cg{new SimulationModel<Bead>};
        cg->readFrom(stream, catalog_);
        return cg;
    }
  
}

