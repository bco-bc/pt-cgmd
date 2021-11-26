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
        ff_ptr_t forceField(std::istream& stream, const spec_catalog_ptr_t& catalog);

        /**
         * Obtains a force field from an input file.
         * @param fileName File name.
         * @param catalog Particle specification catalog.
         * @return Force field.
         */
        ff_ptr_t forceField(const std::string& fileName, const spec_catalog_ptr_t& catalog);

        /**
         * Returns current force field.
         * @return Force field.
         */
        ff_ptr_t forceField();

        /**
         * Returns a default simulation parameters.
         * @return  Simulation parameters. Holds default values for 'simulation.nsteps',
         * 'simulation.nwrite', 'simulation.npairlists', 'simulation.temperature', and
         * 'simulation.timestep'.
         */
        sim_param_ptr_t simulationParameters();

        /**
         * Returns pair lists generator.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Pair lists generator.
         */
        pair_list_gen_ptr_t pairListsGenerator(const bc_ptr_t &bc);

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
                                    const bc_ptr_t &bc);
        
        /**
         * Returns an algorithm (a displacer) to displace a particle system (to change
         * the system's state).
         * @param displacerType Displacer identifier.
         * @return Algorithm.
         */
        displacer_ptr_t displacer(std::string displacerType,
                                  const sim_param_ptr_t& simulationParameters,
                                  const interactor_ptr_t& interactor);
        
        /**
         * Returns periodic boundary conditions.
         * @param box Simulation box.
         * @return Boundary condition.
         */
        bc_ptr_t 
        boundaryCondition(const box_ptr_t& box);
        
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
        forces_ptr_t forces(const bc_ptr_t& bc, const ff_ptr_t& forceField);

        /**
         * Reads a particle system from an input file.
         * @param fileName File name.
         * @param catalog Particle specifications catalog.
         * @param isCoarseGrained Specifies whether the input file holds a coarse-grained particle system.
         * @return Particle system.
         */
        p_system_ptr_t particleSystem(std::string fileName,
                                      const spec_catalog_ptr_t& catalog,
                                      bool isCoarseGrained);

    }
}

#endif /* S_FACTORY_HPP */

