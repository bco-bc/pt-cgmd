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
    class pair_lists_generator;
    class Interactor;
    class PairLists;
    class pair_potential;
    class RF;
    class SF;
    class SC;
    class Forces;
    class displacer;

    // S is protonation state type, e.g. Discrete or Continuous.
    template <typename S>
    class ProtonatableBead;

    // S is protonation state type, e.g. Discrete or Continuous.
    template <typename S>
    class ProtonatableCoarseGrained;

    // Types.

    /*
     * Boundary condition.
     */
    using bc_t = boundary_condition;

    /**
     * Boundary condition pointer typeName.
     */
    using bc_ptr_t = std::shared_ptr<bc_t>;

    /**
     * Force field pointer typeName.
     */
    using ff_ptr_t = std::shared_ptr<ForceField>;

    /**
     * Type for holding parameters, such as values for the time step
     * and the requested number of steps.
     */
    using param_t = param::param_t;
    using param_ptr_t = param::param_ptr_t;

    /**
     * Protonatable coarse grained particle system.
    */
    using prot_cg_sys_t = ProtonatableCoarseGrained<Continuous>;
    // using prot_cg_sys_t = ProtonatableCoarseGrained<Discrete>;

    /**
     * Protonatable coarse grained particle system pointer.
     */
    using prot_cg_sys_ptr_t = std::shared_ptr<prot_cg_sys_t>;

    /**
     * Protonatable bead typeName.
     */
    using prot_bead_t = ProtonatableBead<Continuous>;
    //using prot_bead_ptr_t = ProtonatableCoarseGrainedBead<Discrete>;

    /**
     * Protonatable bead pointer typeName.
     */
    using prot_bead_ptr_t = std::shared_ptr<prot_bead_t>;

    /**
     * Protonatable particle system factory pointer typeName.
     */
    using prot_p_sys_factory = std::shared_ptr<ProtonatableParticleSystemFactory>;

    /**
     * Pair list generator typeName.
     */
    using pair_list_gen_ptr_t = std::shared_ptr<pair_lists_generator>;

    /**
     * Force calculation pointer typeName.
     */
    using forces_ptr_t = std::shared_ptr<Forces>;

    /**
     * Pair potential pointer typeName.
     */
    using pair_potential_ptr_t = std::shared_ptr<pair_potential>;

    /**
     * Interactor pointer typeName.
     */
    using interactor_ptr_t = std::shared_ptr<Interactor>;

    /**
     * displacer pointer typeName.
     */
    using displacer_ptr_t = std::shared_ptr<displacer>;

    /**
     * Reaction field electrostatic interaction pointer typeName.
     */
    using rf_ptr_t = std::shared_ptr<RF>;

    /**
     * Shifted force electrostatic interaction pointer typeName.
     */
    using sf_ptr_t = std::shared_ptr<SF>;

    /**
     * Screened Coulomb electrostatic interaction pointer typeName.
     */
    using sc_ptr_t = std::shared_ptr<SC>;

    namespace units {

        /**
         * Converter between molecular units and MVV_DPD units.
         * @tparam V
         */
        template <typename V>
        class dpd;

        using dpd_ptr_t = std::shared_ptr<dpd<real_t>>;
    }



    // OLD

    class SimulationModelFactory;
    class ProtonTransferPairListGenerator;
    class ProtonTransfer;


    // P is particle typeName, e.g. Atom or Bead.
    template <typename P>
    class SimulationModel;

    /**
     * Coarse grained simulation model.
     */
    //using cg_sim_model_ptr_t = std::shared_ptr<SimulationModel<Bead>>;
    
    //using sim_model_fact_ptr_t = std::shared_ptr<SimulationModelFactory>;
            

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
     * Protonatable beads pair list generator pointer typeName.
     */
    using prot_pair_list_gen_ptr_t = std::shared_ptr<ProtonTransferPairListGenerator>;
    
    /**
     * Proton transfer displacer pointer typeName.
     */
    using pt_displacer_ptr_t = std::shared_ptr<ProtonTransfer>;


}

#endif /* S_TYPES_HPP */

