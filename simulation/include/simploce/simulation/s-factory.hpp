/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 3:39 PM
 */

#ifndef S_FACTORY_HPP
#define S_FACTORY_HPP

#include "s-types.hpp"
#include "simploce/particle/p-factory.hpp"


namespace simploce {
    namespace factory {

        /**
         * Obtains a force field from an input stream.
         * @param stream Input stream.
         * @param catalog Particle specification catalog.
         * @return Force field.
         */
        ff_ptr_t obtainFrom(std::istream& stream, const spec_catalog_ptr_t& catalog);

        /**
         * Returns current force field.
         * @return Force field.
         */
        ff_ptr_t forceField();

        /**
         * Returns a default simulation parameters.
         * @return  Simulation parameters.
         */
        sim_param_ptr_t simulationParameters();

        /**
         * Returns pair lists generator.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Pair lists generator.
         */
        pair_list_gen_ptr_t pairListsGenerator(const box_ptr_t &box,
                                               const bc_ptr_t &bc);

        /**
         * Returns coarse grained force field for harmonic potential 
         * applicable to bonds between beads. 
         * @param catalog Particle specification catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @param fc Force constant harmonic potential.
         * @param Rref Reference distance.
         * @return Force field.
         *
        cg_ff_ptr_t
        harmonicPotentialForceField(const spec_catalog_ptr_t& catalog,
                                    const bc_ptr_t& bc,
                                    const box_ptr_t& box,
                                    real_t fc,
                                    const length_t& Rref);
        
        /**
         * Coarse grained force field for polarizable water.
         * @param catalog Particle specification catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @param protonatable If true, water is protonatable.
         * @return Force field.
         *
        cg_ff_ptr_t
        polarizableWaterForceField(const spec_catalog_ptr_t& catalog,
                                   const bc_ptr_t& bc,
                                   const box_ptr_t& box,
                                   bool protonatable = false);
        
        /**
         * Returns coarse grained force field for formic acid (HCOOH) in water.
         * @param catalog Particle specification catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @param water Coarse grained force field for water.
         * @return Force field.
         *
        cg_ff_ptr_t 
        formicAcidSolutionForceField(const spec_catalog_ptr_t& catalog,
                                     const bc_ptr_t& bc,
                                     const box_ptr_t& box,
                                     const cg_ff_ptr_t& water);
        
        /**
         * Returns coarse grained force field for an NaCl electrolyte solution.
         * @param catalog Particle specification catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @return Force field.
         *
        cg_ff_ptr_t
        electrolyteForceField(const spec_catalog_ptr_t& catalog,
                              const bc_ptr_t& bc,
                              const box_ptr_t& box);
        
        /**
         * Returns coarse grained force field for a LJ fluid.
         * @param catalog Particle specification catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @return  Force field.
         *
        cg_ff_ptr_t
        ljFluidForceField(const spec_catalog_ptr_t& catalog,
                          const bc_ptr_t& bc,
                          const box_ptr_t& box);
         */

        /**
         * Returns factory for creating protonatable particle systems.
         * @return Factory.
         */
        prot_p_sys_factory protonatableParticleSystemFactory(const spec_catalog_ptr_t& catalog);

        /**
         * Returns an interactor.
         * @return Interactor.
         */
        interactor_ptr_t interactor(const sim_param_ptr_t& simulationParameters,
                                    const ff_ptr_t& forceField,
                                    const box_ptr_t &box,
                                    const bc_ptr_t &bc);
        
        /**
         * Returns simulation model factory.
         * @param particle_model_fact_ptr_t Particle model factory.
         * @param catalog Particle specifications catalog.
         * @return Model factory.
         *
        sim_model_fact_ptr_t 
        simulationModelFactory(const spec_catalog_ptr_t& catalog);
        
        /**
         * Returns coarse grained interactor for a harmonic potential force field
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @param fc Force constant harmonic potential.
         * @param Rref Reference distance.
         * @return Interactor.
         *
        cg_interactor_ptr_t
        harmonicPotentialInteractor(const spec_catalog_ptr_t& catalog,
                                    const box_ptr_t& box, 
                                    const bc_ptr_t& bc,
                                    real_t fc,
                                    const length_t& Rref);
        
        
        /**
         * Returns coarse grained interactor.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @param protonatable If true, the interactor will assume parameters for 
         * protonatable waters.
         * @return Interactor.
         *
        cg_interactor_ptr_t 
        polarizableWaterInteractor(const spec_catalog_ptr_t& catalog,
                                   const box_ptr_t& box, 
                                   const bc_ptr_t& bc,
                                   bool protonatable = false);
        
        /**
         * Returns coarse grained interactor.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @return Interactor.
         *
        cg_interactor_ptr_t
        formicAcidSolutionInteractor(const spec_catalog_ptr_t& catalog,
                                     const box_ptr_t& box, 
                                     const bc_ptr_t& bc);
        
        /**
         * Returns coarse grained interactor for electrolyte solution.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @return Interactor.
         *
        cg_interactor_ptr_t
        electrolyteInteractor(const spec_catalog_ptr_t& catalog,
                              const box_ptr_t& box, 
                              const bc_ptr_t& bc);
        
        /**
         * Returns coarse grained interactor for LJ fluid.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @return Interactor.
         *
        cg_interactor_ptr_t
        ljFluidInteractor(const spec_catalog_ptr_t& catalog,
                              const box_ptr_t& box, 
                              const bc_ptr_t& bc);
        
        /**
         * Leap frog algorithm for atomistic particle models.
         * @param interactor Atomistic interactor.
         * @return Algorithm.
         *
        at_displacer_ptr_t 
        leapFrog(at_interactor_ptr_t& interactor);
        
        /**
         * Leap frog algorithm for coarse grained particle models.
         * @param interactor Coarse grained interactor.
         * @return Algorithm.
         *
        cg_displacer_ptr_t 
        leapFrog(cg_interactor_ptr_t& interactor);
        
        /**
         * Velocity Verlet algorithm for atomistic particle models.
         * @param interactor Atomistic interactor.
         * @return Algorithm.
         *
        at_displacer_ptr_t 
        velocityVerlet(at_interactor_ptr_t& interactor);
        
        /**
         * Velocity Verlet algorithm for coarse grained particle models.
         * @param interactor Coarse grained interactor.
         * @return Algorithm.
         *
        cg_displacer_ptr_t 
        velocityVerlet(cg_interactor_ptr_t& interactor);

        /**
         * Langevin Velocity Verlet algorithm for atomistic particle models.
         * @param interactor Atomistic interactor.
         * @return Algorithm.
         *
        at_displacer_ptr_t 
        langevinVelocityVerlet(at_interactor_ptr_t& interactor);
        
        /**
         * Langevin Velocity Verlet algorithm for coarse grained particle models.
         * @param interactor Coarse grained interactor.
         * @return Algorithm.
         *
        cg_displacer_ptr_t 
        langevinVelocityVerlet(cg_interactor_ptr_t& interactor);
        
        /**
         * Returns Langevin Velocity Verlet algorithm with proton transfer for 
         * coarse grained models.
         * @param interactor Coarse grained interactor.
         * @param generator Protonatable bead pair list generator.
         * @param displacer Proton transfer displacer.
         * @return Displacer.
         *
        cg_displacer_ptr_t
        protonTransferlangevinVelocityVerlet(const cg_interactor_ptr_t& interactor,
                                             const prot_pair_list_gen_ptr_t& generator,
                                             const pt_displacer_ptr_t& displacer);
        
        /**
         * Change displacer according to given specification.
         * @param displacerSpec Displacer specification.
         * @param sm Coarse grained simulation model.
         *
        void
        changeDisplacer(const std::string& displacerSpec, cg_sim_model_ptr_t& sm);
        */

        /**
         * Returns periodic boundary conditions.
         * @param box Simulation box.
         * @return Boundary condition.
         */
        bc_ptr_t 
        pbc(const box_ptr_t& box);
        
        /**
         * Returns generator of pairs of protonatable beads possibly involved in proton
         * transfer.
         * @param bc Boundary condition.
         * @return Generator.
         *
        prot_pair_list_gen_ptr_t
        protonTransferPairListGenerator(const bc_ptr_t& bc);
        
        /**
         * Returns proton transfer (PT) displacer with constant rate.
         * @param rate Rate.
         * @param gamma Inverse of time constant of decay.
         * @return PT displacer
         *
        pt_displacer_ptr_t 
        protonTransferDisplacer();
         */

        /**
         * Returns a protonatable coarse-grained particle system
         * @return Particle system.
         */
        prot_cg_sys_ptr_t protonatableCoarseGrained();

        /**
         * Returns force calculator.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @param forceField Force field.
         * @return Force calculator.
         */
        forces_ptr_t forces(const box_ptr_t& box, const bc_ptr_t& bc, const ff_ptr_t& forceField);
        
    }
}

#endif /* FACTORY_HPP */

