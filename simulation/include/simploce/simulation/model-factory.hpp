/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   model-factory.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 29, 2019, 5:22 PM
 */

#ifndef MODEL_FACTORY_HPP
#define MODEL_FACTORY_HPP

#include "cg-displacer.hpp"
#include "sim-model.hpp"
#include "stypes.hpp"
#include "simploce/particle/coarse-grained.hpp"

namespace simploce {
    
    class ModelFactory {
    public:
        
        /**
         * Constructor.
         * @param catalog Particle specifications catalog.
         */
        ModelFactory(const spec_catalog_ptr_t& catalog);
    
        /**
         * Returns a coarse-grained (CG) simulation model of polarizable water 
         * according to Riniker and van Gunsteren, J. Chem. Phys. 134, 084119, 2011. 
         * The model maps 5 waters onto a polarizable coarse-grained bead with two mass points.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param atDensitySI Atomistic density in SI units kg/m^3. Default is 
         * 997.0479 kg/m^3 at 298.15 K (25 C).
         * @param temperature Temperature, default is 298.15 K.
         * @return Coarse grained model.
         * @see <a href="https://en.wikipedia.org/wiki/Density#Water">Water Density</a>
         */
        cg_sim_model_ptr_t createPolarizableWater(const spec_catalog_ptr_t& catalog,
                                                  const box_ptr_t& box,
                                                  const density_t atDensitySI = 997.0479,
                                                  const temperature_t temperature = 298.15);
        
        /**
         * Create a coarse grained simulation of formic acid in water.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param molarity Molarity. Default 0.1 M.
         * @param temperature Temperature. Default is 298.15 K.
         * @return Coarse grained model.
         */
        cg_sim_model_ptr_t createFormicAcidSolution(const spec_catalog_ptr_t& catalog,
                                                    const box_ptr_t& box,
                                                    const molarity_t molarity = 0.1,
                                                    const temperature_t temperature = 298.15);
        
        /**
         * Creates new coarse grained model by reading it from an input stream.
         * @param stream Input stream.
         * @return Coarse grained model.
         */
        cg_sim_model_ptr_t readCoarseGrainedFrom(std::istream& stream);
        
    private:
        
        spec_catalog_ptr_t catalog_;
    };
}


#endif /* MODEL_FACTORY_HPP */

