/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 3, 2019, 12:45 PM
 */

#ifndef PARTICLE_SYSTEM_FACTORY_HPP
#define PARTICLE_SYSTEM_FACTORY_HPP

#include "p-types.hpp"
#include "p-factory.hpp"

namespace simploce {

    /**
     * Create particle model.
     */
    class ParticleSystemFactory {
    public:
        
        /**
         * Constructor
         * @param catalog Particle specifications catalog.
         */
        explicit ParticleSystemFactory(spec_catalog_ptr_t  catalog);

        /**
         * Create a diatomic molecule consisting of two identical atom type (e.g. O).
         * @param distance Distance between atoms.
         * @param spec Particle specification of atoms.
         * @param temperature Temperature (K). Velocities are assigned to the atoms compatible
         * with the given temperature.
         * @return Diatomic molecule.
         */
        p_system_ptr_t diatomic(const dist_t& distance,
                                const spec_ptr_t& spec,
                                const temperature_t& temperature = 298.15);

        /**
         * Returns coarse grained polarizable water model. This model is based on
         * Riniker and van Gunsteren, J. Chem. Phys., 134, 084110, 2011.
         * @param box Box. Default side length is 7.27 nm.
         * @param nLimit Upper limit for the number of CG waters. One CG water consists of two water
         * beads, CW and DP.
         * @param densitySI Requested atomistic density in SI units (kg/m^3).
         * Default is 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param temperature Requested temperature in SI units (K). Default is
         * 298.15 K.
         * @return Polarizable water.
         * @see <a href="https://aip.scitation.org/doi/10.1063/1.3553378">
         * Riniker and van Gunsteren</a>
         */
        p_system_ptr_t polarizableWater(const box_ptr_t& box = factory::box(7.27),
                                        std::size_t nLimit = 2560,
                                        const density_t& densitySI = 997.0479,
                                        const temperature_t& temperature = 298.15);

        /**
         * Returns electrolyte solution of given molarity in a box. This model
         * is loosely based upon the Debye-Huckel theory of electrolytes.
         * Ions are represented by Na+ and Cl- particles, and water is
         * considered as a background continuum. The Debye-Huckel theory is
         * valid for low ionic strength only.
         * @param box Particle box
         * @param molarity Molarity in mol/l (M). Default is 0.1 M. The average
         * physiological salt concentration is 0.15 mM.
         * @param temperature Requested temperature in K. Default is 298.15 K.
         * @return Ionic solution.
         */
        p_system_ptr_t simpleElectrolyte(const box_ptr_t& box = factory::box(6.30),
                                         const molarity_t& molarity = 0.1,
                                         const temperature_t& temperature = 298.15);

        /**
         * Returns argon at given density and temperature. The defaults values are
         * compatible with liquid argon as employed by Rahman in 1964.
         * @param box Simulation box.
         * @param nLimit Maximal number of atoms allowed.
         * @param densitySI Requested atomistic density in SI units (kg/m^3).
         * Default is 1.374 g/cm^3 = 1374.0 kg/m^3 (at a pressure of 1 atm).
         * @param temperature Requested temperature in SI units (K). Default is
         * 94.4 K.
         * @return Atomistic particle model.
         * @see <a href="https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405">Rahman</a>
         */
        p_system_ptr_t argon(const box_ptr_t& box = factory::box(3.47786),
                             std::size_t nLimit = 864,
                             const density_t& densitySI = 1374.0,
                             const temperature_t& temperature = 94.4);

    protected:

        velocity_t generateVelocity_(const temperature_t& temperature);

        /**
         * Returns particle specification catalog.
         * @return Catalog
         */
        spec_catalog_ptr_t& catalog();

    private:
        
        spec_catalog_ptr_t catalog_;
    };
}

#endif /* PARTICLE_SYSTEM_FACTORY_HPP */

