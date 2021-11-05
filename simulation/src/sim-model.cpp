/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 12:07 PM
 */

#include "simploce/simulation/sim-model.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/cg-displacer.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/no-bc.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/protonatable-bead.hpp"
#include "simploce/util/cube.hpp"
#include "simploce/simulation/constant-rate-pt.hpp"
#include <memory>
#include <utility>
#include <stdexcept>
#include <iostream>

namespace simploce {
    
    SimulationModel<Bead>::SimulationModel() :
        cg_{}, displacer_{}, interactor_{}, box_{}, bc_{}
    {        
    }

    SimulationModel<Bead>::SimulationModel(const cg_mod_ptr_t& cg,
                                           const cg_displacer_ptr_t& displacer,                        
                                           const cg_interactor_ptr_t& interactor,
                                           const box_ptr_t& box,
                                           const bc_ptr_t& bc) :
       cg_{cg}, displacer_{displacer}, interactor_{interactor}, box_{box}, bc_{bc}
    {
    }
       
    std::size_t 
    SimulationModel<Bead>::size() const
    {
        return cg_->numberOfParticles();
    }
       
    SimulationData 
    SimulationModel<Bead>::interact(const sim_param_t& param)
    {
        SimulationData data;
        auto result = interactor_->interact(param, cg_);
        data.bepot = result.first;
        data.nbepot = result.second;
        return data;
    }
    
    SimulationData 
    SimulationModel<Bead>::interact(const bead_ptr_t& bead,
                                    const sim_param_t& param)
    {
        SimulationData data;
        auto result = interactor_->interact(bead, param, cg_);
        data.bepot = result.first;
        data.nbepot = result.second;
        return data;
    }
    
    SimulationData 
    SimulationModel<Bead>::displace(const sim_param_t& param)
    { 
        return displacer_->displace(param, cg_); 
    }
        
    void 
    SimulationModel<Bead>::saveState(std::ostream& stream)
    {
        /*
        cg_->doWithAll<void>([&stream] (const std::vector<bead_ptr_t>& beads) {
            for (auto bead: beads) {
                bead->writeState(stream);
            }
        });
        cg_->doWithProtonatableBeads<void>([&stream] (const std::vector<prot_bead_ptr_t>& beads) {
            for (auto const bead : beads) {
                bead->writeState(stream);
            }
        });
         */
    }
    
    void 
    SimulationModel<Bead>::readState(std::istream& stream)
    {
        /*
         cg_->doWithAll<void>([&stream] (const std::vector<bead_ptr_t>& beads) {
            for (auto bead: beads) {
                bead->readState(stream);
            }
        });
        cg_->doWithProtBeads<void>([&stream] (const std::vector<dprot_bead_ptr_t>& discrete,
                                              const std::vector<cprot_bead_ptr_t>& continuous) {
            for (auto d : discrete) {
                d->readState(stream);
            }
            for (auto c : continuous) {
                c->readState(stream);
            }
        });
         */
    }
    
    void 
    SimulationModel<Bead>::readFrom(std::istream& stream,
                                    const spec_catalog_ptr_t& catalog)
    {
        std::string stringBuffer;
        
        // Particles.
        cg_ = CoarseGrained::obtainFrom(stream, catalog);
        
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
            interactor_ = factory::polarizableWaterInteractor(catalog, box_, bc_);
        } else if (stringBuffer.find(conf::ELECTROLYTE, 0) != std::string::npos) {
            interactor_ = factory::electrolyteInteractor(catalog, box_, bc_);
        } else if (stringBuffer.find(conf::LJ_FLUID, 0) != std::string::npos) {
            interactor_ = factory::ljFluidInteractor(catalog, box_, bc_);
        } else {
            throw std::domain_error("Cannot identify force field.");
        }
        
        // Displacer
        std::getline(stream, stringBuffer);
        if ( stringBuffer.find(conf::LANGEVIN_VELOCITY_VERLET, 0) != std::string::npos ) {
            displacer_ = factory::langevinVelocityVerlet(interactor_);
        } else if ( stringBuffer.find(conf::LEAP_FROG, 0) != std::string::npos) {
            displacer_ = factory::leapFrog(interactor_);
        } else {
            throw std::domain_error("Cannot identify displacer.");
        }        
    }
    
    void 
    SimulationModel<Bead>::displacer(const cg_displacer_ptr_t& displacer)
    {
        if ( !displacer ) {
            throw std::domain_error("No displacer specified.");
        }
        displacer_ = displacer;
        std::clog << "Displacer set to " << displacer->id() << std::endl;
    }
    
    cg_interactor_ptr_t 
    SimulationModel<Bead>::interactor() const
    {
       return interactor_; 
    }
    
    cg_displacer_ptr_t 
    SimulationModel<Bead>::displacer() const
    {
        return displacer_;
    }
    
    std::ostream& 
    operator << (std::ostream& stream, 
                 const SimulationModel<Bead>& sm)
    {
        stream << *sm.cg_ << std::endl;
        stream << *sm.box_ << std::endl;
        stream << sm.bc_->id() << std::endl;
        stream << sm.interactor_->id() << std::endl;
        stream << sm.displacer_->id();
        return stream;
    }

}
