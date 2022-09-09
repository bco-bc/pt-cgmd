/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 3, 2019, 12:45 PM
 */

#ifndef PARTICLE_SYSTEM_FACTORY_HPP
#define PARTICLE_SYSTEM_FACTORY_HPP

#include "p-types.hpp"
#include "p-factory.hpp"
#include "simploce/conf/p-conf.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/plane.hpp"

namespace simploce {

    /**
     * Create particle model. Particle position are always initially in the positive quadrant of the
     * Cartesian coordinate system.
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
         * @param temperature Requested temperature (K). Velocities are assigned to the atoms compatible
         * with the given temperature.
         * @return Diatomic molecule.
         */
        p_system_ptr_t diatomic(const dist_t& distance,
                                const spec_ptr_t& spec,
                                const temperature_t& temperature = units::si<real_t>::ROOM_TEMPERATURE);

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
         * is loosely based upon the Debye-Hückel theory of electrolytes.
         * Ions are represented by Na+ and Cl- particles, and water is
         * considered as a background continuum. The Debye-Hückel theory is
         * valid for low ionic strength only.
         * @param box Particle box
         * @param molarity Molarity in mol/l (M). Default is 0.1 M. The average
         * physiological salt concentration is 0.15 mM.
         * @param temperature Requested temperature in K. Default is room temperature.
         * @return Ionic solution.
         */
        p_system_ptr_t simpleElectrolyte(const box_ptr_t& box = factory::box(6.30),
                                         const molarity_t& molarity = 0.1,
                                         const temperature_t& temperature = units::si<real_t>::ROOM_TEMPERATURE);

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

        /**
         * Creates a coarse-grained (that is, mesoscale) polymer solution in water.
         * @param box Box
         * @param chainLength Numbers of monomeric units (beads) in the chain, each of which will be of
         * the same particle specification and represent several real monomers.
         * @param monomericUnitBeadSpecName Specification name for the monomeric units.
         * @param numberOfPolymers Required number of polymers.
         * @param spacing Spacing between polymer segments in a chain.
         * @param numberOfWaters Required number of water beads.
         * @waterBeadSpecName Specification name for water beads.
         * @poram temperature Requested temperature. This assumes a DPD unit system, where temperature is
         * expressed as kT=n, where n is a non-negative integer and k is the Boltzmann constant.
         * @param placeRandom Places particles randomly in the box. By default, a grid is employed to create
         * the polymer solution.
         * @return Particle system.
         */
        p_system_ptr_t polymerSolution(const box_ptr_t& box,
                                       int chainLength,
                                       const std::string& monomericUnitBeadSpecName,
                                       int numberOfPolymers,
                                       const length_t& spacing,
                                       int numberOfWaters,
                                       const std::string& waterBeadSpecName,
                                       const temperature_t& temperature = 1,
                                       bool placeRandom = false);


        /**
         * Creates a coarse-grained (that is, mesoscale) mixture of droplet, polymer and solvent beads
         * (the latter is referred to as water). This follows the setup by Howard et al, Soft Matter, 2019, 15, 3168.
         * @param box Box
         * @param chainLength Numbers of monomeric units (beads) in the chain, each of which will be of
         * the same particle specification and represent several real monomers.
         * @param monomericUnitBeadSpecName Specification name for the monomeric units.
         * @param numberOfPolymers Required number of polymers.
         * @param spacing Spacing between polymer segments in a chain.
         * @param numberOfWaters Required number of water (solvent) beads.
         * @param waterBeadSpecName Specification name for water beads.
         * @param numberOfDropletBeads Number of droplet particles (beads).
         * @param dropletBeadSpecName Droplet bead name.
         * @param temperature Requested temperature. This assumes a DPD unit system, where temperature is
         * expressed as kT=n, where n is a non-negative integer and k is the Boltzmann constant.
         * @return Particle system.
         */
        p_system_ptr_t dropletPolymerSolution(const box_ptr_t& box,
                                              int chainLength,
                                              const std::string& monomericUnitBeadSpecName,
                                              int numberOfPolymers,
                                              const length_t& spacing,
                                              int numberOfWaters,
                                              const std::string& waterBeadSpecName,
                                              int numberOfDropletBeads,
                                              const std::string& dropletBeadSpecName,
                                              const temperature_t& temperature = 1);


        /**
         * Adds a layer of surface/boundary particles.
         * @param particleSystem Particle system to which boundary particles are added.
         * @param spacing Spacing between boundary particles.
         * @param plane Plane identifying surface.  Both corresponding sides of the particle
         * box will be covered with boundary particles.
         */
        void addParticleBoundary(const p_system_ptr_t& particleSystem,
                                 dist_t spacing,
                                 Plane plane,
                                 bool excludeCorner = false);

        /**
         * Generates a atomic particle system from PDB.
         * @param fileName Filename with PDB entry.
         * @param temperature Requested temperature in K.
         * @param excludeWater If true, water in the input file is ignored.
         * @return Particle system.
         */
        p_system_ptr_t fromPDB(std::string& fileName,
                               const temperature_t& temperature = units::si<real_t>::ROOM_TEMPERATURE,
                               bool excludeWater = true);

    protected:

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

