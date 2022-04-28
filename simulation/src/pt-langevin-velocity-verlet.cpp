/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 25, 2019, 2:47 PM
 */

#include "simploce/simulation/pt-langevin-velocity-verlet.hpp"
#include "simploce/simulation/pt-pair-list-generator.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/simulation/pt.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include <stdexcept>

namespace simploce {
    
    static cg_displacer_ptr_t lvv_{};
    
    ProtonTransferLangevinVelocityVerlet::
        ProtonTransferLangevinVelocityVerlet(const cg_interactor_ptr_t& interactor,
                                             const pt_pair_list_gen_ptr_t& generator,
                                             const pt_displacer_ptr_t& displacer) : 
    interactor_{interactor}, generator_{generator}, displacer_{displacer}
    {        
        if ( !interactor ) {
            throw std::domain_error(
                "PT LangevinVelocityVerlet: Missing interactor."
            );
        }
        if ( !generator ) {
            throw std::domain_error(
                "PT LangevinVelocityVerlet: Missing protonatable bead pair list generator."
            );
        }
        if ( !displacer_ ) {
            throw std::domain_error(
                "PT LangevinVelocityVerlet: Missing displacer for proton transfer."
            );
        }
        lvv_ = factory::langevinVelocityVerlet(interactor_);
    }
    
    SimulationData 
    ProtonTransferLangevinVelocityVerlet::displace(const sim_param_t& param, 
                                                   const cg_mod_ptr_t& cg) const
    {
        using prot_pair_list_t = ProtonTransferPairListGenerator::prot_pair_list_t;
        
        static std::size_t counter = 0;
        static prot_pair_list_t pairlist{};
        
        counter += 1;
        
        std::size_t npairlists = param.get<std::size_t>("npairlists");

        /*
        if ( counter % npairlists == 0 ) {
            pairlist = generator_->generate(cg);
        }
        
        // Transfer protons to update mass and charge values.
        if ( !pairlist.empty() ) {
            cg->doWithProtonatableBeads<void>([this, param] (const std::vector<prot_bead_ptr_t>& beads) {
                this->displacer_->transfer(param, beads, pairlist);
            });
        }
         */
        
        // Update positions and velocities.
        SimulationData data = lvv_->displace(param, cg);
        data.numberOfProtonTransferPairs = pairlist.size();
                
        return data;
    }
    
    std::string 
    ProtonTransferLangevinVelocityVerlet::id() const
    {
        return conf::PT_LANGEVIN_VELOCITY_VERLET;
    }
}

