/*
 * Author: André H. Juffer, Biocenter Oulu, University of Finland.
 *
 * Created on October 28, 2021.
 */

#ifndef SIMULATION_PROTONATABLE_PARTICLE_SYSTEM_FACTORY_HPP
#define SIMULATION_PROTONATABLE_PARTICLE_SYSTEM_FACTORY_HPP

#include "simploce/types/s-types.hpp"
#include "s-factory.hpp"
#include "simploce/particle/particle-system-factory.hpp"


namespace simploce {

    class ProtonatableParticleSystemFactory: public ParticleSystemFactory {
    public:

        explicit ProtonatableParticleSystemFactory(const spec_catalog_ptr_t& catalog);

        /**
         * Returns a coarse grained protonatable polarizable particle model for water. This model
         * extends the model by Riniker and van Gunsteren, J. Chem. Phys., 134, 084110, 2011.
         * @param box Box.
         * @param nLimit Upper limit for the number of CG water
         * @param densitySI Requested atomistic density in SI units (kg/m^3).
         * Default is 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param temperature Requested temperature in SI units (K). Default is
         * 298.15 K.s.
         * @return Protonatable coarse grained particle model.
         * @see <a href="https://aip.scitation.org/doi/10.1063/1.3553378">
         * Riniker and van Gunsteren</a>
         */
        prot_cg_sys_ptr_t
        protonatablePolarizableWater(const box_ptr_t& box = factory::box(7.27),
                                     std::size_t nLimit = 2560,
                                     const density_t& densitySI = 997.0479,
                                     const temperature_t& temperature = 298.15);

        /**
         * Returns coarse grained protonatable coarse grained model for a formic acid solution,
         * consisting of protonatable beads. Water is also polarizable.
         * @param box Box.
         * @param nLimitWater Maximum number of CG waters.
         * @param densitySI Requested atomistic water density in SI units (kg/m^3).
         * Default is 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param molarity Requested molarity of formic acid. Default is 0.1 M.
         * @param temperature Requested temperature in SI units (K). Default is
         * 298.15 K.
         * @return Protonatable coarse grained particle model.
         * @see #protonatablePolarizableWater
         */
        prot_cg_sys_ptr_t formicAcid(const box_ptr_t& box = factory::box(7.27),
                                     std::size_t nLimitWater = 2560,
                                     const density_t& densitySI = 997.0479,
                                     const molarity_t& molarity = 0.1,
                                     const temperature_t& temperature = 298.15);


    };
}

#endif //SIMULATION_PROTONATABLE_PARTICLE_SYSTEM_FACTORY_HPP
