/*
 * Factory methods.
 * File:   p-factory.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 4:56 PM
 */

#ifndef P_FACTORY_HPP
#define P_FACTORY_HPP

#include "p-types.hpp"
#include <string>
#include <iostream>

namespace simploce {
    namespace factory {
        
        /**
         * Creates a -cubic- box.
         * @param side Side length.
         * @return Box. Always new instance.
         */
        box_ptr_t box(const length_t& side);

        /**
         * Obtains a particle specification catalog from a file.
         * @param fileName Input file name.
         * @return Particle specification catalog.
         */
        spec_catalog_ptr_t particleSpecCatalog(const std::string& fileName);
        
        /**
         * Obtains a particle specification catalog from an input stream.
         * @param stream Input stream.
         * @return Particle specification catalog.
         */
        spec_catalog_ptr_t particleSpecCatalog(std::istream& stream);
        
        /**
         * Returns particle model factory.
         * @param catalog Particle specification catalog.
         * @return Particle model factory.
         */
        p_system_fact_ptr_t particleModelFactory(const spec_catalog_ptr_t& catalog);
        
        /**
         * Returns an empty atomistic particle model.
         * @return Particle model.
         */
        at_sys_ptr_t atomistic();
        
        /**
         * Returns an empty coarse grained particle model.
         * @return Particle model.
         */
        cg_sys_ptr_t coarseGrained();
        
    }
}


#endif /* P_FACTORY_HPP */

