/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 3:41 PM
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/potentials/forces.hpp"
#include "simploce/simulation/distance-pair-list-generator.hpp"
#include "simploce/types/s-types.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/mc.hpp"
#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/1d-pbc.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/util/file.hpp"
#include <memory>
#include <stdexcept>

namespace simploce {
    namespace factory {

        // Force field.
        static ff_ptr_t forceField_{};

        // Forces calculator.
        static forces_ptr_t forces_{};

        // Particle pair list generator.
        static pair_list_gen_ptr_t pairListsGenerator_{};

        // Simulation parameters.
        static param_ptr_t simulationParameters_{};
        
        // Protonatable particle system factory.
        static prot_p_sys_factory protonatableParticleModelFactory_{};

        // Simulation model factory.
        //static sim_model_fact_ptr_t simModelFactory_{};

        // Boundary conditions.
        static bc_ptr_t bc_{};
        
        // Interactor.
        static interactor_ptr_t interactor_{};

        //static prot_pair_list_gen_ptr_t ptPairlisGen_{};
        
        //static pt_displacer_ptr_t constantRate_{};


        ff_ptr_t
        forceField(std::istream& stream, const spec_catalog_ptr_t& catalog) {
            if ( !forceField_ ) {
                forceField_ = ForceField::obtainFrom(stream, catalog);
            }
            return forceField_;
        }

        ff_ptr_t forceField(const std::string& fileName, const spec_catalog_ptr_t& catalog) {
            if ( !forceField_) {
                std::ifstream stream;
                util::open_input_file(stream, fileName);
                forceField_ = factory::forceField(stream, catalog);
                stream.close();
            }
            return forceField_;
        }

        ff_ptr_t forceField() {
            if ( !forceField_ ) {
                throw std::domain_error("No force field defined");
            }
            return forceField_;
        }

        param_ptr_t
        simulationParameters() {
            if ( !simulationParameters_ ) {
                simulationParameters_ = std::make_shared<simploce::param::param_t>();
                simulationParameters_->put("simulation.nsteps", 1000);
                simulationParameters_->put("simulation.nwrite", "10");
                simulationParameters_->put("simulation.npairlists", 10);
                simulationParameters_->put("simulation.temperature", 298.15);  // K.
                simulationParameters_->put("simulation.timestep", 0.001);      // 1 fs.
                simulationParameters_->put("simulation.gamma", 1.0); // ps^-1.
                simulationParameters_->put("simulation.include-external", 0);  // False
                simulationParameters_->put("forces.nb.cutoff", 2.6);
                simulationParameters_->put("displacer.mc.range", 0.1);
            }
            return simulationParameters_;
        }

        pair_list_gen_ptr_t pairListsGenerator(const dist_t& cutoff, const bc_ptr_t &bc) {
            if ( !pairListsGenerator_ ) {
                pairListsGenerator_ = std::make_shared<DistancePairListGenerator>(cutoff, bc);
            }
            return pairListsGenerator_;
        }

        prot_p_sys_factory protonatableParticleSystemFactory(const spec_catalog_ptr_t& catalog) {
            if ( !protonatableParticleModelFactory_ ) {
                protonatableParticleModelFactory_ =
                        std::make_shared<ProtonatableParticleSystemFactory>(catalog);
            }
            return protonatableParticleModelFactory_;
        }


        displacer_ptr_t displacer(const std::string& displacerType,
                                  const param_ptr_t& simulationParameters,
                                  const interactor_ptr_t& interactor) {
            static util::Logger logger("simploce::factory::displacer");
            if (displacerType == conf::LEAP_FROG ) {
                logger.debug("Creating LeapFrog displacer.");
                return std::make_shared<LeapFrog>(simulationParameters, interactor);
            } else if ( displacerType == conf::MONTE_CARLO ) {
                 logger.debug("Creating Monte Carlo displacer.");
                return std::make_shared<MonteCarlo>(simulationParameters, interactor);
            } else if ( displacerType == conf::LANGEVIN_VELOCITY_VERLET ) {
                logger.debug("Creating Langevin Velocity Verlet displacer.");
                return std::make_shared<LangevinVelocityVerlet>(simulationParameters, interactor);
            } else if ( displacerType == conf::VELOCITY_VERLET ) {
                logger.debug("Creating Velocity Verlet displacer.");
                return std::make_shared<VelocityVerlet>(simulationParameters, interactor);
            } else {
                util::logAndThrow(logger, displacerType + ": No such displacer.");
                return nullptr;
            }
        }

        /*
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
        boundaryCondition(const box_ptr_t& box)
        {
            if ( !bc_ ) {
                bc_ = std::make_shared<PeriodicBoundaryCondition>(box);
            }
            return bc_;
        }

        bc_ptr_t
        oneDimensionBoundaryCondition(const box_ptr_t& box, const Direction& direction)
        {
            if ( !bc_ ) {
                bc_ = std::make_shared<OneD_PeriodicBoundaryCondition>(box, direction);
            }
            return bc_;
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
        interactor_ptr_t interactor(const param_ptr_t& simulationParameters,
                                    const ff_ptr_t& forceField,
                                    const bc_ptr_t &bc) {
            if ( !interactor_ ) {
                dist_t cutoff = simulationParameters->get<real_t>("forces.nb.cutoff");
                auto pairListGenerator = factory::pairListsGenerator(cutoff, bc);
                auto forces = factory::forces(simulationParameters, bc, forceField);
                interactor_ = std::make_shared<Interactor>(simulationParameters,
                                                           pairListGenerator,
                                                           forces);
            }
            return interactor_;
        }


        prot_cg_sys_ptr_t protonatableCoarseGrained() {
            return std::make_shared<prot_cg_sys_t>();
        }

        forces_ptr_t forces(const param_ptr_t& simulationParam, const bc_ptr_t& bc, const ff_ptr_t& forceField) {
            if ( !forces_ ) {
                forces_ = std::make_shared<Forces>(simulationParam, bc, forceField);
            }
            return forces_;
        }

        p_system_ptr_t particleSystem(const std::string& fileName,
                                      const spec_catalog_ptr_t& catalog,
                                      bool isCoarseGrained) {
            std::ifstream stream;
            util::open_input_file(stream, fileName);
            p_system_ptr_t particleSystem;
            if ( isCoarseGrained ) {
                particleSystem = prot_cg_sys_t::obtainFrom(stream, catalog);
            } else {
                auto ps = Atomistic::obtainFrom(stream, catalog);
                particleSystem = ps;
            }
            stream.close();
            return particleSystem;
        }
    }
}
