/*
 * Author: André H. Juffer, Biocenter Oulu, University of Finland.
 *
 * Created on October 28, 2021.
 */

#ifndef SIMULATION_PROTONATABLE_PARTICLE_MODEL_FACTORY_HPP
#define SIMULATION_PROTONATABLE_PARTICLE_MODEL_FACTORY_HPP

#include "s-types.hpp"
#include "simploce/particle/particle-model-factory.hpp"


namespace simploce {

    class ProtonatableParticleModelFactory: public ParticleModelFactory {
    public:

        explicit ProtonatableParticleModelFactory(const spec_catalog_ptr_t& catalog);

        /**
         * Returns a coarse grained protonatable polarizable particle model for water. This model
         * extends the model by Riniker and van Gunsteren, J. Chem. Phys., 134, 084110, 2011.
         * @param box Box.
         * @param densitySI Requested atomistic density in SI units (kg/m^3).
         * Default is 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param temperature Requested temperature in SI units (K). Default is
         * 298.15 K.
         * @param nLimit Upper limit for the number of CG waters.
         * @return Protonatable coarse grained particle model.
         * @see <a href="https://aip.scitation.org/doi/10.1063/1.3553378">
         * Riniker and van Gunsteren</a>
         */
        prot_cg_mod
        protonatablePolarizableWater(const box_ptr_t& box,
                                     const density_t& densitySI = 997.0479,
                                     const temperature_t& temperature = 298.15,
                                     std::size_t nLimit = 1000000);

        /**
         * Returns coarse grained protonatable coarse grained model for a formic acid solution,
         * consisting of protonatable beads.
         * @param box Box.
         * @param densitySI Requested atomistic water density in SI units (kg/m^3).
         * Default is 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param molarity Requested molarity of formic acid. Default is 0.1 M.
         * @param temperature Requested temperature in SI units (K). Default is
         * 298.15 K.
         * @return Protonatable coarse grained particle model.
         */
        prot_cg_mod formicAcid(const box_ptr_t& box,
                               const density_t& densitySI = 997.0479,
                               const molarity_t& molarity = 0.1,
                               const temperature_t& temperature = 298.15);


    };
}

#endif //SIMULATION_PROTONATABLE_PARTICLE_MODEL_FACTORY_HPP
