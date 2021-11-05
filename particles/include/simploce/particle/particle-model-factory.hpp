/*
 * To change this license header, choose License Headers in Project Properties.
 * To protonate this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   particle-model-factory.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 3, 2019, 12:45 PM
 */

#ifndef PARTICLE_MODEL_FACTORY_HPP
#define PARTICLE_MODEL_FACTORY_HPP

#include "atomistic.hpp"
#include "coarse-grained.hpp"
#include "protonatable-coarse-grained.hpp"
#include "p-types.hpp"

namespace simploce {

    /**
     * Create particle model.
     */
    class ParticleModelFactory {
    public:
        
        /**
         * Constructor
         * @param catalog Particle specifications catalog.
         */
        explicit ParticleModelFactory(const spec_catalog_ptr_t& catalog);

        /**
         * Create a diatomic molecule consisting of two identical atom type (e.g. O).
         * @param distance Distance between atoms.
         * @param spec Particle specification of atoms.
         * @param temperature Temperature (K). Velocities are assigned to the atoms compatible
         * with the given temperature.
         * @return Diatomic molecule.
         */
        Atomistic diatomic(const distance_t& distance,
                           const spec_ptr_t& spec,
                           const temperature_t& temperature = 298.15);

        /**
         * Returns coarse grained polarizable water model. This model is based on
         * Riniker and van Gunsteren, J. Chem. Phys., 134, 084110, 2011.
         * @param box Box.
         * @param densitySI Requested atomistic density in SI units (kg/m^3).
         * Default is 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param temperature Requested temperature in SI units (K). Default is 
         * 298.15 K.
         * @param nLimit Upper limit for the number of CG waters.
         * @return Coarse grained particle model.
         * @see <a href="https://aip.scitation.org/doi/10.1063/1.3553378">
         * Riniker and van Gunsteren</a>
         */
        CoarseGrained polarizableWater(const box_ptr_t& box,
                                       const density_t& densitySI = 997.0479,
                                       const temperature_t& temperature = 298.15,
                                       std::size_t nLimit = 1000000);
        
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
         * @return Coarse grained particle model.
         */
        Atomistic simpleElectrolyte(const box_ptr_t& box,
                                    const molarity_t& molarity = 0.1,
                                    const temperature_t& temperature = 298.15);
        
        /**
         * Returns argon at given density and temperature. Defaults values are compatible with liquid argon.
         * @param box Simulation box.
         * @param densitySI Requested atomistic density in SI units (kg/m^3).
         * Default is 1.374 g/cm^3 = 1374.0 kg/m^3 (at a pressure of 1 atm).
         * @param temperature Requested temperature in SI units (K). Default is 
         * 94.4 K.
         * @return Atomistic particle model.
         * @see <a href="https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405">Argon at Wikipedia</a>
         */
        Atomistic argon(const box_ptr_t& box,
                        const density_t& densitySI = 1374.0,
                        const temperature_t& temperature = 94.4);

    protected:

        /**
         * Returns particle specification catalog.
         * @return Catalog
         */
        spec_catalog_ptr_t catalog();

        /**
         * Assigns/sets the velocity of the given atom according to the given temperature.
         * @param atom Atom.
         * @param temperature Temperature.
         */
        void assignVelocity(atom_ptr_t& atom, const temperature_t& temperature);

        /**
         * Assigns/sets the velocity of the given bead according to the given temperature.
         * @param bead bead.
         * @param temperature Temperature.
         */
        void assignVelocity(bead_ptr_t& bead, const temperature_t& temperature);
        
    private:
        
        spec_catalog_ptr_t catalog_;
    };
}

#endif /* PARTICLE_MODEL_FACTORY_HPP */

