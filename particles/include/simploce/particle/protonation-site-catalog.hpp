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
 * File:   protonation-site-catalog.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 12, 2019, 5:29 PM
 */

#ifndef PROTONATION_SITE_CATALOG_HPP
#define PROTONATION_SITE_CATALOG_HPP

#include "atomistic.hpp"
#include "coarse-grained.hpp"
#include "ptypes.hpp"
#include <iostream>

namespace simploce {
    
    /**
     * Holds specifications of protonation sites.
     */
    class ProtonationSiteCatalog {
    public:
        
        /**
         * Looks up protonation sites composed of atoms.
         * @param at Atomistic particle model.
         * @return Protonation sites.
         */
        std::vector<atom_prot_site_ptr_t> lookup(Atomistic& at) const;
        
        /**
         * Returns all protonation site names.
         * @return Names.
         */
        std::vector<std::string> names() const;
        
        /**
         * Returns catalog.
         * @param stream Input stream from which protonation sites are obtained.
         * @return Catalog.
         */
        static prot_site_catalog_ptr_t create(std::istream& stream);
        
    private:
        
        ProtonationSiteCatalog();                
        
    };
        
    /**
     * Writes protonation site catalog to output stream.
     * @param stream Output stream.
     * @param catalog Protonation site catalog.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const ProtonationSiteCatalog& catalog);
}

#endif /* PROTONATION_SITE_CATALOG_HPP */

