/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 8, 2019, 2:23 PM
 */

#ifndef PARTICLE_CATALOG_HPP
#define PARTICLE_CATALOG_HPP

#include "p-types.hpp"
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
         * Creates a particle specification catalog from given particle specifications.
         * @param specs Particle specifications.
         * @return Catalog
         */
        static spec_catalog_ptr_t
        create(std::map<std::string, spec_ptr_t>& specs);

        /**
         * Inquires whether the catalog has a particle specification of given name.
         * @param name Particle specification name
         * @return True if it has, otherwise false.
         */
        bool hasSpecification(const std::string& name) const;
        
        /**
         * Returns a specification.
         * @param name Specification name. Must not contain spaces.
         * @return Specification.
         * @throw Exception if specification cannot be found.
         */
        spec_ptr_t lookup(const std::string& name) const;

        /**
         * Returns a specification.
         * @param name Element name. Must not contain spaces.
         * @return Specification.
         */
        spec_ptr_t lookupByElementName(const std::string& name) const;
        
        /**
         * Returns specification for molecular water.
         * @return Specification.
         *
         */
        spec_ptr_t molecularWater() const;
        
        /**
         * Returns specification for hydrogen.
         * @return Specification.
         */
        spec_ptr_t H() const;
        
        /**
         * Returns specification for carbon
         * @return Specification.
         */
        spec_ptr_t C() const;
        
        /**
         * Returns specification for oxygen.
         * @return Specification. 
         */
        spec_ptr_t O() const;
        
        /**
         * Returns specification for nitrogen.
         * @return Specification.
         */
        spec_ptr_t N() const;

        /**
         * Returns specification for a static boundary particle.
         * @return Specification, name SBP.
         */
        spec_ptr_t staticBP() const;

        /**
         * Creates a catalog by reading specifications from a given input stream.
         * @param stream Input stream.
         * @return Specifications catalog.
         */
        static spec_catalog_ptr_t obtainFrom(std::istream& stream);
        
        /**
         * Writes this catalog to an output stream.
         * @param stream Output stream.
         */
        void write(std::ostream& stream) const;
        
    private:
        
        using map_specs_t = std::map<std::string, spec_ptr_t>;
                
        explicit ParticleSpecCatalog(map_specs_t specs);
        
        map_specs_t specs_{};
    };
    
    std::ostream&
    operator << (std::ostream& stream,
                 const ParticleSpecCatalog& catalog);
}

#endif /* PARTICLE_CATALOG_HPP */

