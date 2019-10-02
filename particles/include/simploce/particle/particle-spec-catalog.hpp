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
 * File:   particle-catalog.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 8, 2019, 2:23 PM
 */

#ifndef PARTICLE_CATALOG_HPP
#define PARTICLE_CATALOG_HPP

#include "ptypes.hpp"
#include <iostream>
#include <string>
#include <map>

namespace simploce {
    
    /**
     * Holds particle specifications.
     */
    class ParticleSpecCatalog {
    public:
        
        /**
         * Returns a specification.
         * @param name Specification name.
         * @return Specification.
         */
        spec_ptr_t lookup(const std::string& name) const;
        
        /**
         * Returns specification for molecular water.
         * @return Specification.
         */
        spec_ptr_t molecularWater() const;
        
        /**
         * Hydrogen.
         * @return Specification.
         */
        spec_ptr_t H() const;
        
        /**
         * Carbon
         * @return Specification.
         */
        spec_ptr_t C() const;
        
        /**
         * Oxygen.
         * @return Specification. 
         */
        spec_ptr_t O() const;
        
        /**
         * Nitrogen.
         * @return Specification.
         */
        spec_ptr_t N() const;        
        
        /**
         * Creates a catalog by reading specifications from a given input stream.
         * @param stream Input stream.
         * @return Specifications catalog.
         */
        static spec_catalog_ptr_t create(std::istream& stream);
        
        /**
         * Writes this catalog to an output stream.
         * @param stream Output stream.
         */
        void write(std::ostream& stream) const;
        
    private:
        
        using map_specs_t = std::map<std::string, spec_ptr_t>;
                
        ParticleSpecCatalog(const map_specs_t& specs);
        
        map_specs_t specs_{};
    };
    
    std::ostream& operator << (std::ostream& stream, 
                               const ParticleSpecCatalog& catalog);
}

#endif /* PARTICLE_CATALOG_HPP */

