/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 3:41 PM
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/potentials/forces.hpp"
// #include "simploce/simulation/distance-pair-list-generator.hpp"
#include "simploce/simulation/cell-pair-list-generator.hpp"
#include "simploce/types/s-types.hpp"
#include "simploce/units/units-dpd.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/mc.hpp"
#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/mvv-dpd.hpp"
#include "simploce/simulation/s1-dpd.hpp"
#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/pbc-2d.hpp"
#include "simploce/simulation/pbc-1d-bb.hpp"
#include "simploce/simulation/pbc-1d-sr.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/logger.hpp"
#include <memory>
#include <stdexcept>

namespace simploce {
    namespace factory {

        // Force field.
        static ff_ptr_t forceField_{};

        // Forces calculator.
        static forces_ptr_t forces_{};

        // Particle pair list generator.
        static pair_list_gen_ptr_t pairListGenerator_{};

        // Simulation parameters.
        static param_ptr_t param_{};
        
        // Protonatable particle system factory.
        static prot_p_sys_factory protonatableParticleModelFactory_{};

        // Simulation model factory.
        //static sim_model_fact_ptr_t simModelFactory_{};

        // Boundary conditions.
        static bc_ptr_t bc_{};
        
        // Interactor.
        static interactor_ptr_t interactor_{};

        //static prot_pair_list_gen_ptr_t ptPairlistGen_{};
        
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

        ff_ptr_t
        forceField() {
            if ( !forceField_ ) {
                throw std::domain_error("No force field defined");
            }
            return forceField_;
        }

        param_ptr_t
        simulationParameters() {
            if ( !param_ ) {
                param_ = std::make_shared<simploce::param::param_t>();
                param_->put("simulation.nsteps", 1000);
                param_->put("simulation.nwrite", "10");
                param_->put("simulation.npairlists", 10);
                param_->put("simulation.temperature", 298.15);  // K.
                param_->put("simulation.timestep", 0.001);      // 1 fs.
                param_->put("simulation.forces.include-external", 0);  // False
                param_->put("simulation.forces.cutoffSR", 1.8);
                param_->put("simulation.forces.cutoffLR", 2.6);
                param_->put("simulation.displacer.mc.range", 0.1);
                param_->put("simulation.displacer.lvv.gamma", 1.0); // ps^-1
                param_->put("simulation.displacer.dpd.gamma", 3.0);
            }
            return param_;
        }

        pair_list_gen_ptr_t
        pairListGenerator(const param_ptr_t& param, const bc_ptr_t &bc) {
            if ( !pairListGenerator_ ) {
                //pairListGenerator_ = std::make_shared<DistancePairListGenerator>(param, bc);
                pairListGenerator_ = std::make_shared<CellPairListGenerator>(param, bc);
            }
            return pairListGenerator_;
        }

        prot_p_sys_factory
        protonatableParticleSystemFactory(const spec_catalog_ptr_t& catalog) {
            if ( !protonatableParticleModelFactory_ ) {
                protonatableParticleModelFactory_ =
                        std::make_shared<ProtonatableParticleSystemFactory>(catalog);
            }
            return protonatableParticleModelFactory_;
        }

        units::dpd_ptr_t
        dpdUnits(const mass_t& mass, const length_t& length, const energy_t& energy) {
            return std::make_shared<units::dpd<real_t>>(mass, length, energy);
        }

        displacer_ptr_t
        displacer(const std::string& displacerType,
                                  const param_ptr_t& param,
                                  const interactor_ptr_t& interactor,
                                  const bc_ptr_t& bc,
                                  const units::dpd_ptr_t& dpdUnits) {
            static util::Logger logger("simploce::factory::displacer");
            if (displacerType == conf::LEAP_FROG ) {
                logger.debug("Creating LeapFrog displacer.");
                return std::make_shared<LeapFrog>(param, interactor);
            } else if ( displacerType == conf::MONTE_CARLO ) {
                 logger.debug("Creating Monte Carlo displacer.");
                return std::make_shared<MonteCarlo>(param, interactor);
            } else if ( displacerType == conf::LANGEVIN_VELOCITY_VERLET ) {
                logger.debug("Creating Langevin Velocity Verlet displacer.");
                return std::make_shared<LangevinVelocityVerlet>(param, interactor);
            } else if ( displacerType == conf::VELOCITY_VERLET ) {
                logger.debug("Creating Velocity Verlet displacer.");
                return std::make_shared<VelocityVerlet>(param, interactor);
            } else if ( displacerType == conf::DPD) {
                logger.debug("Creating Modified Velocity Verlet DPD displacer.");
                return std::make_shared<MVV_DPD>(param_, interactor, bc, dpdUnits);
            } else if ( displacerType == conf::S1_DPD) {
                logger.debug("Creating S1 DPD displacer.");
                return std::make_shared<S1_DPD>(param_, interactor, bc, dpdUnits);
            } else {
                util::logAndThrow(logger, displacerType + ": No such displacer.");
                return nullptr;
            }
        }

        bc_ptr_t 
        pbc(const box_ptr_t& box)
        {
            if ( !bc_ ) {
                bc_ = std::make_shared<PBC>(box);
            }
            return bc_;
        }

        bc_ptr_t
        pbc_2d(const box_ptr_t& box,
               const Direction& d1,
               const Direction& d2) {
            if ( !bc_ ) {
                bc_ = std::make_shared<PBC_2D>(box, d1, d2);
            }
            return bc_;
        }

        bc_ptr_t
        pbc1dBB(const box_ptr_t& box,
                 const Direction& direction)
        {
            if ( !bc_ ) {
                bc_ = std::make_shared<PBC_1D_BB>(box, direction);
            }
            return bc_;
        }

        bc_ptr_t
        pbc1dSR(const box_ptr_t& box,
                const Direction& direction) {
            if (!bc_) {
                bc_ = std::make_shared<PBC_1D_SR>(box, direction);
            }
            return bc_;
        }

        interactor_ptr_t interactor(const param_ptr_t& param,
                                    const ff_ptr_t& forceField,
                                    const bc_ptr_t &bc) {
            if ( !interactor_ ) {
                auto pairListGenerator = factory::pairListGenerator(param, bc);
                auto forces = factory::forces(param, bc, forceField);
                interactor_ = std::make_shared<Interactor>(param,
                                                           pairListGenerator,
                                                           forces);
            }
            return interactor_;
        }


        prot_cg_sys_ptr_t protonatableCoarseGrained() {
            //return std::make_shared<prot_cg_sys_t>();
            throw std::domain_error("Not implemented");
        }

        forces_ptr_t
        forces(const param_ptr_t& param, const bc_ptr_t& bc, const ff_ptr_t& forceField) {
            if ( !forces_ ) {
                forces_ = std::make_shared<Forces>(param, bc, forceField);
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
                particleSystem = CoarseGrained::parseIt(stream, catalog);
            } else {
                particleSystem = Atomistic::parseIt(stream, catalog);
            }
            stream.close();
            return particleSystem;
        }
    }
}
