/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 3:39 PM
 */

#ifndef S_FACTORY_HPP
#define S_FACTORY_HPP

#include "simploce/types/s-types.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/util/direction.hpp"


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
         * @return Simulation parameters. Holds default values for various parameters.
         */
        param_ptr_t simulationParameters();

        /**
         * Returns pair lists generator.
         * @param SImulation parameters.
         * @param bc Boundary condition.
         * @return Pair lists generator.
         */
        pair_list_gen_ptr_t pairListsGenerator(const param_ptr_t& param , const bc_ptr_t &bc);

        /**
         * Returns factory for creating protonatable particle systems.
         * @return Factory.
         */
        prot_p_sys_factory protonatableParticleSystemFactory(const spec_catalog_ptr_t& catalog);

        /**
         * Returns an interactor.
         * @param param Simulation parameters.
         * @param forceField Force field.
         * @param bc Boundary conditions.
         * @return Interactor.
         */
        interactor_ptr_t interactor(const param_ptr_t& param,
                                    const ff_ptr_t& forceField,
                                    const bc_ptr_t &bc);

        /**
         * Returns converter between molecular units (MU) and MVV_DPD units, based on characteristic values for
         * mass, length, and energy.
         * @param mass Mass in MU.
         * @param length Length in MU.
         * @param energy Energy in MU.
         * @return Converter
         */
        units::dpd_ptr_t dpdUnits(const mass_t& mass, const length_t& length, const energy_t& energy);
        
        /**
         * Returns an algorithm (a displacer) to displace a particle system (to change
         * the system's state).
         * @param displacerType displacer identifier.
         * @param param Parameters for conducting the requested simulation
         * (e.g., number of steps, etc.).
         * @param interactor Interactor.
         * @param bc Boundary condition.
         * @param dpdUnits Units convertor for MVV_DPD simulations. Must be present for MVV_DPD simulations.
         * @return Algorithm.
         */
        displacer_ptr_t displacer(const std::string& displacerType,
                                  const param_ptr_t& param,
                                  const interactor_ptr_t& interactor,
                                  const bc_ptr_t& bc,
                                  const units::dpd_ptr_t& dpdUnits = nullptr);
        
        /**
         * Returns periodic boundary conditions.
         * @param box Simulation box.
         * @return Boundary condition.
         */
        bc_ptr_t 
        boundaryCondition(const box_ptr_t& box);


        /**
         * Returns 1 dimensional periodic boundary conditions.
         * @param box Simulation box.
         * @param direction Apply PBC in this direction only.
         * @return Boundary condition.
         */
        bc_ptr_t
        oneDimensionBoundaryCondition(const box_ptr_t& box, const Direction& direction);
        
        /**
         * Returns generator of pairs of protonatable beads possibly involved in proton
         * transfer.
         * @param bc Boundary condition.
         * @return Generator.
         *
        prot_pair_list_gen_ptr_t
        protonTransferPairListGenerator(const bc_ptr_t& bc);
        
        *
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
         * @param param Simulation parameters.
         * @param bc Boundary condition.
         * @param forceField Force field.
         * @return Force calculator.
         */
        forces_ptr_t forces(const param_ptr_t& param, const bc_ptr_t& bc, const ff_ptr_t& forceField);

        /**
         * Reads a particle system from an input file.
         * @param fileName File name.
         * @param catalog Particle specifications catalog.
         * @param isCoarseGrained Specifies whether the input file holds a coarse-grained particle system.
         * @return Particle system.
         */
        p_system_ptr_t particleSystem(const std::string& fileName,
                                      const spec_catalog_ptr_t& catalog,
                                      bool isCoarseGrained);

    }
}

#endif /* S_FACTORY_HPP */

