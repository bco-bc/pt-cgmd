/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:39 PM
 */

#ifndef S_TYPES_HPP
#define S_TYPES_HPP

#include "simploce/particle/p-types.hpp"
#include "simploce/util/param.hpp"
#include <tuple>

namespace simploce {
    
    // Forward declarations.
    class boundary_condition;
    class ForceField;
    class ProtonatableParticleSystemFactory;
    class Continuous;
    class Discrete;

    // P is particle type, e.g. Atom or Bead.
    template <typename P>
    class pair_lists_generator;

    // S is protonation state type, e.g. Discrete or Continuous.
    template <typename S>
    class ProtonatableBead;

    // S is protonation state type, e.g. Discrete or Continuous.
    template <typename S>
    class ProtonatableCoarseGrained;

    // P is particle type, e.g. Atom or Bead.
    template <typename P>
    class Interactor;

    template <typename P>
    class PairLists;

    // Types.

    /*
     * Boundary condition.
     */
    using bc_t = boundary_condition;

    /**
     * Boundary condition pointer type.
     */
    using bc_ptr_t = std::shared_ptr<bc_t>;

    /**
     * Force field pointer type.
     */
    using ff_ptr_t = std::shared_ptr<ForceField>;

    /**
     * Type for holding simulation parameters, such as values for the time step
     * and the requested number of steps.
     */
    using sim_param_t = param::param_t;

    /**
     * Simulation parameters pointer type.
     */
    using sim_param_ptr_t = std::shared_ptr<sim_param_t>;

    /**
     * Atom pair lists generator pointer type.
     */
    using atom_pair_lists_gen_ptr_t = std::shared_ptr<pair_lists_generator<Atom>>;

    /**
     * Bead pair lists generator pointer type.
     */
    using bead_pair_lists_gen_ptr_t = std::shared_ptr<pair_lists_generator<Bead>>;

    /**
     * Protonatable coarse grained particle system.
    */
    using prot_cg_sys = ProtonatableCoarseGrained<Continuous>;
    // using prot_cg_sys = ProtonatableCoarseGrained<Discrete>;

    /**
     * Protonatable coarse grained particle system pointer.
     */
    using prot_cg_mod_ptr_t = std::shared_ptr<prot_cg_sys>;

    /**
     * Protonatable bead type.
     */
    using prot_bead_t = ProtonatableBead<Continuous>;
    //using prot_bead_ptr_t = ProtonatableCoarseGrainedBead<Discrete>;

    /**
     * Protonatable bead pointer type.
     */
    using prot_bead_ptr_t = std::shared_ptr<prot_bead_t>;

    /**
     * Protonatable particle system factory pointer type.
     */
    using prot_p_sys_factory = std::shared_ptr<ProtonatableParticleSystemFactory>;





    // OLD

    //class CoarseGrainedForceField;
    class CoarseGrainedDisplacer;
    class AtomisticDisplacer;
    class AtomisticForceField;
    class SimulationModelFactory;
    class ProtonTransferPairListGenerator;
    class ProtonTransfer;


    // P is particle type, e.g. Atom or Bead.
    template <typename P>
    class SimulationModel;









    using cg_displacer_ptr_t = std::shared_ptr<CoarseGrainedDisplacer>;
    
    using at_displacer_ptr_t = std::shared_ptr<AtomisticDisplacer>;
    
    using at_ff_ptr_t = std::shared_ptr<AtomisticForceField>;
    
    //using cg_ff_ptr_t = std::shared_ptr<CoarseGrainedForceField>;
    
    //using at_interactor_t = Interactor<Atom>;
        
    //using at_interactor_ptr_t = std::shared_ptr<at_interactor_t>;

    //using cg_interactor_t = Interactor<Bead>;
        
    //using cg_interactor_ptr_t = std::shared_ptr<cg_interactor_t>;
    

    /**
     * Atom pair list generator type.
     */
    using at_ppair_list_gen_t = pair_lists_generator<Atom>;
    

    /**
     * Coarse grained simulation model type.
     */
    using cg_sim_model_t = SimulationModel<Bead>;
    
    /**
     * Coarse grained simulation model.
     */
    using cg_sim_model_ptr_t = std::shared_ptr<SimulationModel<Bead>>;
    
    using sim_model_fact_ptr_t = std::shared_ptr<SimulationModelFactory>;
            

    /**
     * Defines a range of something in a given indexed collection to forms pairs between 
     * elements of the collection.
     * tuple_t t;
     * t.get<0> is the begin index of an indexed collection.
     * t.get<1> is the end index.
     * t.get<2> is the number of pairs.
     */
    //using range_t = std::tuple<std::size_t, std::size_t, std::size_t>;
    
    /**
     * Protonatable beads pair list generator pointer type.
     */
    using prot_pair_list_gen_ptr_t = std::shared_ptr<ProtonTransferPairListGenerator>;
    
    /**
     * Proton transfer displacer pointer type.
     */
    using pt_displacer_ptr_t = std::shared_ptr<ProtonTransfer>;
    
}

#endif /* S_TYPES_HPP */

