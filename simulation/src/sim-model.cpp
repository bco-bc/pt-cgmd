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
 * File:   sim-model.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 12:07 PM
 */

#include "simploce/simulation/sim-model.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/simulation/cg-displacer.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/no-bc.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/protonatable-bead.hpp"
#include "simploce/util/cube.hpp"
#include <memory>
#include <utility>
#include <stdexcept>

namespace simploce {
    
    SimulationModel<Bead>::SimulationModel() :
        cg_{}, displacer_{}, interactor_{}, box_{}, bc_{}
    {        
    }

    SimulationModel<Bead>::SimulationModel(const cg_ptr_t& cg,
                                           const cg_displacer_ptr_t& displacer,                        
                                           const cg_interactor_ptr_t interactor,
                                           const box_ptr_t& box,
                                           const bc_ptr_t& bc) :
       cg_{cg}, displacer_{displacer}, interactor_{interactor}, box_{box}, bc_{bc}
    {
    }
       
    std::size_t SimulationModel<Bead>::size() const
    {
        return cg_->numberOfParticles();
    }
       
    SimulationData SimulationModel<Bead>::interact(const sim_param_t& param)
    {
        SimulationData data;
        data.epot = interactor_->interact(param, cg_);
        return data;
    }
    
    SimulationData SimulationModel<Bead>::displace(const sim_param_t& param)
    { 
        return displacer_->displace(param, cg_); 
    }
        
    void SimulationModel<Bead>::saveState(std::ostream& stream)
    {
        cg_->doWithAll<void>([&stream] (const std::vector<bead_ptr_t>& beads) {
            for (auto bead: beads) {
                bead->writeState(stream);
            }
        });
        cg_->doWithProtBeads<void>([&stream] (const std::vector<prot_bead_ptr_t>& pbeads) {
            for (auto const pbead : pbeads) {
                pbead->writeState(stream);
            }
        });
    }
    
    void SimulationModel<Bead>::readState(std::istream& stream)
    {
         cg_->doWithAll<void>([&stream] (const std::vector<bead_ptr_t>& beads) {
            for (auto bead: beads) {
                bead->readState(stream);
            }
        });
        cg_->doWithProtBeads<void>([&stream] (const std::vector<prot_bead_ptr_t>& pbeads) {
            for (auto pbead : pbeads) {
                pbead->readState(stream);
            }
        });       
    }
    
    void SimulationModel<Bead>::readFrom(std::istream& stream,
                                         const spec_catalog_ptr_t& catalog)
    {
        std::string stringBuffer;
        
        // Particles.
        cg_ = CoarseGrained::readFrom(stream, catalog);
        
        // Simulation box.
        real_t boxSize;
        stream >> boxSize;
        box_ = std::make_shared<box_t>(boxSize);
        std::getline(stream, stringBuffer);  // Read EOL.
        
        // Boundary condition.
        std::getline(stream, stringBuffer);        
        std::size_t n = stringBuffer.find(conf::PBC, 0);
        if (n != std::string::npos ) {
            bc_ = factory::pbc(box_);
        } else {
            std::clog << "Warning: No boundary conditions assumed." << std::endl;
            bc_ = std::make_shared<NoBoundaryCondition>();
        }
        
        // Interactor.
        std::getline(stream, stringBuffer);
        if ( stringBuffer.find(conf::POLARIZABLE_WATER, 0) != std::string::npos ) {
            interactor_ = factory::interactorCoarseGrainedPolarizableWater(catalog, 
                                                                           box_, 
                                                                           bc_);
        } else {
            throw std::domain_error("Cannot identify force field.");
        }
        
        // Displacer
        std::getline(stream, stringBuffer);
        if ( stringBuffer.find(conf::LANGEVIN_VELOCITY_VERLET, 0) != std::string::npos ) {
            displacer_ = factory::langevinVelocityVerlet(interactor_);
        }
        
    }
    
    std::ostream& operator << (std::ostream& stream, const SimulationModel<Bead>& sm)
    {
        stream << *sm.cg_ << std::endl;
        stream << *sm.box_ << std::endl;
        stream << sm.bc_->id() << std::endl;
        stream << sm.interactor_->id() << std::endl;
        stream << sm.displacer_->id();
        return stream;
    }

}
