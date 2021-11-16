/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 3:41 PM
 */

#include "simploce/simulation/s-factory.hpp"
// #include "simploce/simulation/cg-pol-water.hpp"
// #include "simploce/simulation/acid-base-solution.hpp"
// #include "simploce/simulation/cg-electrolyte.hpp"
// #include "simploce/simulation/cg-lj-fluid.hpp"
// #include "simploce/simulation/pair-list-generator.hpp"
// #include "simploce/simulation/cell-lists.hpp"
// #include "simploce/simulation/distance-lists.hpp"
#include "simploce/simulation/distance-pair-list-generator.hpp"
#include "simploce/simulation/s-types.hpp"
// #include "simploce/simulation/leap-frog.hpp"
// #include "simploce/simulation/velocity-verlet.hpp"
// #include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/interactor.hpp"
// #include "simploce/simulation/sim-model-factory.hpp"
// #include "simploce/simulation/pt-pair-list-generator.hpp"
// #include "simploce/simulation/pt-langevin-velocity-verlet.hpp"
#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/simulation/pbc.hpp"
// #include "simploce/simulation/constant-rate-pt.hpp"
//#include "simploce/simulation/cg-hp.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/p-factory.hpp"
#include <memory>
#include <stdexcept>

namespace simploce {
    namespace factory {
        /*
        using at_leap_frog_t = LeapFrog<Atomistic>;
        using cg_leap_frog_t = LeapFrog<CoarseGrained>;
        using at_vv_t = VelocityVerlet<Atomistic>;
        using cg_vv_t = VelocityVerlet<CoarseGrained>;
        using at_lvv_t = LangevinVelocityVerlet<Atomistic>;
        using cg_lvv_t = LangevinVelocityVerlet<CoarseGrained>;
         */
        
        // Force fields.
        static ff_ptr_t forceField_{};                   // Force field.

        //static cg_ff_ptr_t cgPolWaterFF_{};             // Polarizable water.
        //static cg_ff_ptr_t cgFormicAcidSolutionFF_{};   // Formic acid in water.
        //static cg_ff_ptr_t cgElectrolyteFF_{};          // NaCl solution.
        //static cg_ff_ptr_t cgLJFluid_{};                // LJ fluid.
        //static cg_ff_ptr_t cgHPFF_;                      // Harmonic potentials.
        
        // Pair lists generators.
        static bead_pair_lists_gen_ptr_t beadPairlistsGen_{};
        static atom_pair_lists_gen_ptr_t atomPairListsGen_{};

        // Simulation parameters.
        static sim_param_ptr_t simulationParameters_{};
        
        // Leap Frog
        static at_displacer_ptr_t atLeapFrog_{};
        static cg_displacer_ptr_t cgLeapFrog_{};
        
        // Velocity Verlet
        static at_displacer_ptr_t atVV_{};
        static cg_displacer_ptr_t cgVV_{};
        
        // Langevin Velocity Verlet
        static at_displacer_ptr_t atLVV_{};
        static cg_displacer_ptr_t cgLVV_{};
        
        // Langevin Proton Transfer Langevin Velocity Verlet
        static cg_displacer_ptr_t cgPTLVV_{};

        // Protonatable particle model factory.
        static prot_p_sys_factory protonatableParticleModelFactory_{};

        // Simulation model factory.
        static sim_model_fact_ptr_t simModelFactory_{};

        // Boundary conditions.
        static bc_ptr_t pbc_{};
        
        // Interactors.
        static std::shared_ptr<Interactor<Bead>> cgPolWaterInteractor_{};
        static std::shared_ptr<Interactor<Bead>> cgFormicAcidSolutionInteractor_{};
        static std::shared_ptr<Interactor<Bead>> cgElectrolyteInteractor_{};
        static std::shared_ptr<Interactor<Bead>> cgLJFluidInteractor_{};
        static std::shared_ptr<Interactor<Bead>> cgHPInteractor_;
        
        static prot_pair_list_gen_ptr_t ptPairlisGen_{};
        
        static pt_displacer_ptr_t constantRate_{};

        ff_ptr_t
        obtainFrom(std::istream& stream, const spec_catalog_ptr_t& catalog) {
            if ( !forceField_ ) {
                forceField_ = ForceField::obtainFrom(stream, catalog);
            }
            return forceField_;
        }

        ff_ptr_t forceField() {
            if ( !forceField_ ) {
                throw std::domain_error("No force field defined");
            }
            return forceField_;
        }

        sim_param_ptr_t
        simulationParameters() {
            if ( !simulationParameters_ ) {
                simulationParameters_ = std::make_shared<simploce::param::param_t>();
                simulationParameters_->put("simulation.nsteps", 1000);
                simulationParameters_->put("simulation.nwrite", "10");
                simulationParameters_->put("simulation.npairlists", 10);
                simulationParameters_->put("simulation.temperature", 298.15);
            }
            return simulationParameters_;
        }

        bead_pair_lists_gen_ptr_t pairListsGeneratorForBeads(const box_ptr_t &box,
                                                             const bc_ptr_t &bc) {
            if ( !beadPairlistsGen_ ) {
                beadPairlistsGen_ = std::make_shared<DistancePairListGenerator<Bead>>(box, bc);
            }
            return beadPairlistsGen_;
        }

        atom_pair_lists_gen_ptr_t pairListsGeneratorForAtoms(const box_ptr_t &box,
                                                             const bc_ptr_t &bc) {
            if ( !atomPairListsGen_ ) {
                atomPairListsGen_ = std::make_shared<DistancePairListGenerator<Atom>>(box, bc);
            }
            return atomPairListsGen_;
        }

        /*
        cg_ff_ptr_t
        harmonicPotentialForceField(const spec_catalog_ptr_t& catalog,
                                    const bc_ptr_t& bc,
                                    const box_ptr_t& box,
                                    real_t fc,
                                    const length_t& Rref)
        {
            if ( !cgHPFF_ ) {
                cgHPFF_ = 
                    std::make_shared<HarmonicPotential>(catalog, bc, box, fc, Rref);
            }
            return cgHPFF_;
        }
        
                
        cg_ff_ptr_t 
        polarizableWaterForceField(const spec_catalog_ptr_t& catalog,
                                   const bc_ptr_t& bc,
                                   const box_ptr_t& box,
                                   bool protonatable)
        {
            if ( !cgPolWaterFF_ ) {
                cgPolWaterFF_ = 
                    std::make_shared<CoarseGrainedPolarizableWater>(catalog, 
                                                                    bc,
                                                                    box,
                                                                    protonatable);
            }
            return cgPolWaterFF_;
        }
        
        cg_ff_ptr_t 
        formicAcidSolutionForceField(const spec_catalog_ptr_t& catalog,
                                     const bc_ptr_t& bc,
                                     const box_ptr_t& box,
                                     const cg_ff_ptr_t& water)
        {
            if ( !cgFormicAcidSolutionFF_ ) {
                cgFormicAcidSolutionFF_ = 
                    std::make_shared<AcidBaseSolution>(catalog, bc, box, water);
            }
            return cgFormicAcidSolutionFF_;
        }
        
        cg_ff_ptr_t
        electrolyteForceField(const spec_catalog_ptr_t& catalog,
                              const bc_ptr_t& bc,
                              const box_ptr_t& box)
        {
            if ( !cgElectrolyteFF_ ) {
                cgElectrolyteFF_ = 
                    std::make_shared<CoarseGrainedElectrolyte>(catalog, bc, box);
            }
            return cgElectrolyteFF_;
        }
        
        cg_ff_ptr_t
        ljFluidForceField(const spec_catalog_ptr_t& catalog,
                          const bc_ptr_t& bc,
                          const box_ptr_t& box)
        {
            if ( !cgLJFluid_) {
                cgLJFluid_ = 
                    std::make_shared<CoarseGrainedLJFluid>(catalog, bc, box);
            }
            return cgLJFluid_;
        }
        */

        prot_p_sys_factory protonatableParticleModelFactory(const spec_catalog_ptr_t& catalog) {
            if ( !protonatableParticleModelFactory_ ) {
                protonatableParticleModelFactory_ =
                        std::make_shared<ProtonatableParticleSystemFactory>(catalog);
            }
            return protonatableParticleModelFactory_;
        }

        /*
        sim_model_fact_ptr_t 
        simulationModelFactory(const spec_catalog_ptr_t& catalog)
        {
            particle_model_fact_ptr_t particleModelFactory = 
                factory::particleModelFactory(catalog);
            if ( !simModelFactory_ ) {
                simModelFactory_ =
                    std::make_shared<SimulationModelFactory>(particleModelFactory, 
                                                             catalog);
            }
            return simModelFactory_;
        }
        
        cg_interactor_ptr_t
        harmonicPotentialInteractor(const spec_catalog_ptr_t& catalog,
                                    const box_ptr_t& box, 
                                    const bc_ptr_t& bc,
                                    real_t fc,
                                    const length_t& Rref)
        {
            if ( !cgHPInteractor_ ) {
                bead_pair_lists_gen_ptr_t generator =
                    factory::pairListsGeneratorForBeads(box, bc);
                cg_ff_ptr_t forcefield = 
                    factory::harmonicPotentialForceField(catalog, bc, box, fc, Rref);
                cgHPInteractor_ =
                    std::make_shared<Interactor<Bead>>(forcefield, generator);
            }
            return cgHPInteractor_;
        }

        
        cg_interactor_ptr_t 
        polarizableWaterInteractor(const spec_catalog_ptr_t& catalog,
                                   const box_ptr_t& box, 
                                   const bc_ptr_t& bc,
                                   bool protonatable)        
        {
            if ( !cgPolWaterInteractor_ ) {
                bead_pair_lists_gen_ptr_t generator =
                    factory::pairListsGeneratorForBeads(box, bc);
                cg_ff_ptr_t forcefield = 
                    polarizableWaterForceField(catalog, bc, box, protonatable);
                cgPolWaterInteractor_ = 
                    std::make_shared<Interactor<Bead>>(forcefield, generator);
            }
            return cgPolWaterInteractor_;
        }
        
        cg_interactor_ptr_t
        formicAcidSolutionInteractor(const spec_catalog_ptr_t& catalog,
                                     const box_ptr_t& box, 
                                     const bc_ptr_t& bc)
        {
            if ( !cgFormicAcidSolutionInteractor_) {
                bead_pair_lists_gen_ptr_t generator =
                    factory::pairListsGeneratorForBeads(box, bc);
                cg_ff_ptr_t water = 
                    polarizableWaterForceField(catalog, bc, box, true);
                cg_ff_ptr_t forcefield = 
                    formicAcidSolutionForceField(catalog, bc, box, water);
                cgFormicAcidSolutionInteractor_ = 
                    std::make_shared<Interactor<Bead>>(forcefield, generator);
            }
            return cgFormicAcidSolutionInteractor_;
        }
        
        cg_interactor_ptr_t
        electrolyteInteractor(const spec_catalog_ptr_t& catalog,
                              const box_ptr_t& box, 
                              const bc_ptr_t& bc)
        {
            if ( !cgElectrolyteInteractor_) {
                bead_pair_lists_gen_ptr_t generator =
                    factory::pairListsGeneratorForBeads(box, bc);
                cg_ff_ptr_t forcefield =
                    factory::electrolyteForceField(catalog, bc, box);
                cgElectrolyteInteractor_ = 
                    std::make_shared<Interactor<Bead>>(forcefield, generator);
            }
            return cgElectrolyteInteractor_;
        }
        
        cg_interactor_ptr_t
        ljFluidInteractor(const spec_catalog_ptr_t& catalog,
                              const box_ptr_t& box, 
                              const bc_ptr_t& bc)        
        {
            if ( !cgLJFluidInteractor_ ) {
                bead_pair_lists_gen_ptr_t generator =
                    factory::pairListsGeneratorForBeads(box, bc);
                cg_ff_ptr_t forcefield =
                    factory::ljFluidForceField(catalog, bc, box);
                cgLJFluidInteractor_ =
                    std::make_shared<Interactor<Bead>>(forcefield, generator);
            }
            return cgLJFluidInteractor_;
        }
        
        at_displacer_ptr_t 
        leapFrog(at_interactor_ptr_t& interactor)
        {
            if ( !atLeapFrog_ ) {
                atLeapFrog_ = std::make_shared<at_leap_frog_t>(interactor);
            }
            return atLeapFrog_;
        }
        
        cg_displacer_ptr_t 
        leapFrog(cg_interactor_ptr_t& interactor)
        {
            if ( !cgLeapFrog_ ) {
                cgLeapFrog_ = std::make_shared<cg_leap_frog_t>(interactor);
            }
            return cgLeapFrog_;
            
        }
        
        at_displacer_ptr_t 
        velocityVerlet(at_interactor_ptr_t& interactor)
        {
            if ( !atVV_ ) {
                atVV_ = std::make_shared<at_vv_t>(interactor);
            }
            return atVV_;            
        }
        
        cg_displacer_ptr_t 
        velocityVerlet(cg_interactor_ptr_t& interactor)
        {
            if ( !cgVV_ ) {
                cgVV_ = std::make_shared<cg_vv_t>(interactor);
            }
            return cgVV_;                        
        }
        
        at_displacer_ptr_t 
        langevinVelocityVerlet(at_interactor_ptr_t& interactor)
        {
            if ( !atLVV_ ) {
                atVV_ = std::make_shared<at_lvv_t>(interactor);
            }
            return atLVV_;                        
        }
        
        cg_displacer_ptr_t 
        langevinVelocityVerlet(cg_interactor_ptr_t& interactor)
        {
            if ( !cgLVV_ ) {
                cgLVV_ = std::make_shared<cg_lvv_t>(interactor);
            }
            return cgLVV_;                                    
        }
        
        cg_displacer_ptr_t
        protonTransferlangevinVelocityVerlet(const cg_interactor_ptr_t& interactor,
                                             const prot_pair_list_gen_ptr_t& generator,
                                             const pt_displacer_ptr_t& displacer)
        {
            if ( !cgPTLVV_ ) {
                cgPTLVV_ = 
                    std::make_shared<ProtonTransferLangevinVelocityVerlet>(interactor,
                                                                           generator,
                                                                           displacer);
            }
            return cgPTLVV_;
        }   
        
        void 
        changeDisplacer(const std::string& displacerId, cg_sim_model_ptr_t& sm)
        {
            if ( sm->displacer()->id() != displacerId) {                            
                auto interactor = sm->interactor();
                cg_displacer_ptr_t displacer;
                if ( displacerId == conf::LEAP_FROG ) {
                    displacer = factory::leapFrog(interactor);                
                } else if ( displacerId == conf::VELOCITY_VERLET ) {
                    displacer = factory::velocityVerlet(interactor);
                } else if ( displacerId == conf::LANGEVIN_VELOCITY_VERLET ) {
                    displacer = factory::langevinVelocityVerlet(interactor);                
                } else if ( displacerId == conf::PT_LANGEVIN_VELOCITY_VERLET ) {
                    auto bc = sm->boundaryCondition();
                    auto ptGenerator = factory::protonTransferPairListGenerator(bc);
                    auto ptDisplacer = factory::protonTransferDisplacer();
                    displacer = factory::protonTransferlangevinVelocityVerlet(interactor, 
                                                                              ptGenerator, 
                                                                              ptDisplacer);
                } else {
                    throw std::domain_error(displacerId + ": No such displacer.");
                }
                sm->displacer(displacer);
            }
        }
        */

        bc_ptr_t 
        pbc(const box_ptr_t& box)
        {
            if ( !pbc_ ) {
                pbc_ = std::make_shared<PeriodicBoundaryCondition>(box);
            }
            return pbc_;
        }

        /*
        prot_pair_list_gen_ptr_t
        protonTransferPairListGenerator(const bc_ptr_t& bc)
        {
            if ( !ptPairlisGen_) {
                ptPairlisGen_ = 
                    std::make_shared<ProtonTransferPairListGenerator>(bc);
            }
            return ptPairlisGen_;
        }
        
        pt_displacer_ptr_t 
        protonTransferDisplacer()
        {
            if ( !constantRate_ ) {
                constantRate_ = 
                    std::make_shared<ConstantRateProtonTransfer>();
            }
            return constantRate_;
        }
        */
    }
}
