/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 29, 2019, 5:22 PM
 */

#ifndef SIM_MODEL_FACTORY_HPP
#define SIM_MODEL_FACTORY_HPP

#include "cg-displacer.hpp"
#include "sim-model.hpp"
#include "s-types.hpp"
#include "simploce/particle/coarse-grained.hpp"

namespace simploce {
    
    class SimulationModelFactory {
    public:
        
        /**
         * Constructor.
         * @param particleModelFactory Particle model factory.
         * @param catalog Particle specifications catalog.
         */
        SimulationModelFactory(const particle_model_fact_ptr_t& particleModelFactory,
                               const spec_catalog_ptr_t& catalog);
        
        /**
         * Returns coarse-grained (CG) simulation model of consisting of a 
         * single particle group containing one bond whose beads undergo 
         * harmonic motion. Particles are placed along the z-axis.
         * @param R0 Initial distance between particles. Default value is 0.5 nm.
         * @param Rref Reference distance between beads. Default value is 0.4 nm.
         * @param fc Force constant for harmonic potential. Default value is 
         * 100 kJ/(mol nm^2).
         * @return Coarse grained simulation model.
         */
        cg_sim_model_ptr_t harmonic(length_t R0 = 0.5, 
                                    length_t Rref = 0.4, 
                                    real_t fc = 100.0);
                
    
        /**
         * Returns a coarse-grained (CG) simulation model of polarizable water 
         * according to Riniker and van Gunsteren, J. Chem. Phys. 134, 084119, 2011. 
         * The model maps 5 waters onto a polarizable coarse-grained bead with two mass points.
         * @param box Simulation box.
         * @param atDensitySI Atomistic water density in SI units kg/m^3. Default is 
         * 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param temperature Temperature, default is 298.15 K.
         * @param nlimit Upper limit for the number of CG waters.
         * @return Coarse grained simulation model.
         * @see <a href="https://en.wikipedia.org/wiki/Density#Water">Water Density</a>
         */
        cg_sim_model_ptr_t polarizableWater(const box_ptr_t& box,
                                            const density_t atDensitySI = 997.0479,
                                            const temperature_t temperature = 298.15,
                                            std::size_t nlimit = 1000000);
        
        /**
         * Create a coarse grained simulation of formic acid in water.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param atDensitySI Atomistic water density in SI units kg/m^3. Default is 
         * 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param molarity Molarity. Default 0.1 M.
         * @param temperature Temperature. Default is 298.15 K.
         * @param protonatable If true, both water and HCOOH are protonatable, in which
         * case all HCOOH beads will be protonated initially, If false, no 
         * (de)protonation will occur and all HCOOH beads will be deprotonated.
         * @return Coarse grained simulation model.
         */
        cg_sim_model_ptr_t formicAcidSolution(const box_ptr_t& box,
                                              bool protonatable = true,
                                              const molarity_t molarity = 0.1,
                                              const density_t atDensitySI = 997.0479,
                                              const temperature_t temperature = 298.15);
        
        /**
         * Returns electrolyte solution of given molarity in a box. This model 
         * is loosely based upon the Debye-Huckel theory of electrolytes.
         * Ions are represented by Na+ and Cl- particles, while water is 
         * considered as a background continuum. The Debye-Huckel theory is 
         * valid for low ionic strength only. The force field 
         * consists of simple Lennard-Jones and screened Coulomb interactions 
         * between ions. PBC is applied.
         * @param box Simulation box.
         * @param molarity Molarity in mol/l (M). Default is 1.0 M. Note 
         * that physiological salt concentration is about 0.15 mM.
         * @param temperature Requested temperature in K. Default is 298.15 K.
         * @return Coarse grained simulation model.
         */
        cg_sim_model_ptr_t electrolyte(const box_ptr_t& box,
                                       molarity_t molarity = 1.0,
                                       temperature_t temperature = 298.15);
        
        
        /**
         * Returns LJ fluid of given in a simulation box. The force field 
         * consists of simple Lennard-Jones interactions between LJ 
         * beads. PBC is applied.
         * @param box Simulation box.
         * @param densitySI Density in SI units kg/m^3. Default is 
         * 997.0479 kg/m^3 at 298.15 K (25 C) (= atomistic water density).
         * @param temperature Requested temperature in K. Default is 298.15 K.
         * @return Coarse grained simulation model.
         */
        cg_sim_model_ptr_t ljFluid(const box_ptr_t& box,
                                   const density_t densitySI = 997.0479,
                                   const temperature_t temperature = 298.15);
        /**
         * Creates new coarse grained model by reading it from an input stream.
         * @param stream Input stream.
         * @return Coarse grained model.
         */
        cg_sim_model_ptr_t readCoarseGrainedFrom(std::istream& stream);
        
    private:
        
        particle_model_fact_ptr_t particleModelFactory_;
        spec_catalog_ptr_t catalog_;
    };
}


#endif /* MODEL_FACTORY_HPP */

