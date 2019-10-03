/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
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

#include "polarizable-water.hpp"
#include "ptypes.hpp"

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
        ParticleModelFactory(const spec_catalog_ptr_t& catalog);
        
        /**
         * Returns coarse grained polarizable water model.
         * @param box Box.
         * @param atDensitySI Requested atomistic density in SI units (kg/m^3). 
         * Default is 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param temperature Requested temperature in SI units (K). Default is 
         * 298.15 K.
         * @return Coarse grained particle model.
         */
        cg_pol_water_ptr_t polarizableWater(const box_ptr_t& box,
                                            const density_t atDensitySI = 997.0479,
                                            const temperature_t temperature = 298.15) const;
        
        /**
         * Returns coarse grained model for a formic acid solution.
         * @param box Box.
         * @param atDensitySI Requested atomistic water density in SI units (kg/m^3). 
         * Default is 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param molarity Requested molarity of formic acid. Default is 0.1 M.
         * @param temperature Requested temperature in SI units (K). Default is 
         * 298.15 K.
         * @return Coarse grained particle model.
         */
        cg_ptr_t formicAcidSolution(const box_ptr_t& box,
                                    const density_t atDensitySI = 997.0479,
                                    const molarity_t molarity = 0.1,
                                    const temperature_t temperature = 298.15) const;
        
    private:
        
        spec_catalog_ptr_t catalog_;
    };
}

#endif /* PARTICLE_MODEL_FACTORY_HPP */

