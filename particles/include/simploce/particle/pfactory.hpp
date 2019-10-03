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
 * File:   pfactory.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 4:56 PM
 */

#ifndef PFACTORY_HPP
#define PFACTORY_HPP

#include "ptypes.hpp"
#include <string>
#include <iostream>

namespace simploce {
    namespace factory {
        
        /**
         * Returns cubic box.
         * @param side Side length.
         * @return Box. Always new instance.
         */
        box_ptr_t cube(const length_t& side);        
        
        /**
         * Reads a specification catalog from a file.
         * @param fileName Input file name.
         * @return Particle specification catalog.
         */
        spec_catalog_ptr_t particleSpecCatalog(const std::string& fileName);
        
        /**
         * Reads a specification catalog from an input stream.
         * @param stream Input stream.
         * @return Particle specification catalog.
         */
        spec_catalog_ptr_t particleSpecCatalog(std::istream& stream);
        
        /**
         * Returns particle model factory.
         * @param catalog Particle specification catalog.
         * @return Particle model factory.
         */
        particle_model_fact_ptr_t particleModelFactory(const spec_catalog_ptr_t& catalog);
        
        /**
         * Returns an empty atomistic particle model.
         * @return Particle model.
         */
        at_ptr_t atomistic();
        
        /**
         * Returns an empty coarse grained particle model.
         * @return Particle model.
         */
        cg_ptr_t coarseGrained();
        
    }
}


#endif /* PFACTORY_HPP */

