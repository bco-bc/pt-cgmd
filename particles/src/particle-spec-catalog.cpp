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
 * File:   particle-catalog.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 8, 2019, 2:29 PM
 */

#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/pconf.hpp"
#include <boost/algorithm/string.hpp>
#include <stdexcept>
#include <sstream>
#include <fstream>

namespace simploce {
    
    ParticleSpecCatalog::ParticleSpecCatalog(const map_specs_t& specs) : 
       specs_{specs}
    {
    }
    
    particle_spec_ptr_t ParticleSpecCatalog::lookup(const std::string& particleName) const
    {
        using iter_t = map_specs_t::const_iterator;
        
        iter_t iter = specs_.find(particleName);
        if ( iter == specs_.end() ) {
            throw std::domain_error(
                particleName + ": No particle specification available for this particle."
            );
        }
        return (*iter).second;
    }
    
    particle_spec_ptr_t ParticleSpecCatalog::molecularWater() const
    {
        return this->lookup("mH2O");
    }
    
    spec_catalog_ptr_t ParticleSpecCatalog::create(std::istream& stream)
    {
        using pair_t = std::pair<std::string, particle_spec_ptr_t>;
        map_specs_t specs{};
        
        const size_t bufferSize = conf::NAME_WIDTH;
        char charBuffer[bufferSize];
        std::string stringBuffer;
        
        real_t charge, mass, radius, pKa;
        
        // Skip first line.
        std::getline(stream, stringBuffer);
        
        // Read first specification type name ("atom", "bead", "prot_bead").
        stream.read(charBuffer, bufferSize);
        while ( stream.good() ) {
            std::string tname = std::string(charBuffer, bufferSize);            
            boost::trim(tname);
            stream.read(charBuffer, bufferSize);
            std::string name = std::string(charBuffer, bufferSize); 
            boost::trim(name);
            stream >> mass >> charge >> radius;            
            if ( tname == "atom") {
                atom_spec_ptr_t spec = 
                    ParticleSpec::createForAtom(name, charge, mass, radius);
                pair_t p = std::make_pair(name, spec);
                specs.insert(p);
            } else if (tname == "bead") {
                bead_spec_ptr_t spec = 
                    ParticleSpec::createForBead(name, charge, mass, radius);
                pair_t p = std::make_pair(name, spec);
                specs.insert(p);                
            } else if (tname == "prot_bead" ) {
                stream >> pKa;
                prot_bead_spec_ptr_t spec = 
                    ParticleSpec::createForProtonatableBead(
                        name, charge, mass, radius, pKa
                    );
                pair_t p = std::make_pair(name, spec);
                specs.insert(p);                
            } else {
                throw std::domain_error(tname + ": Unknown particle specification type.");
            }
            
            stringBuffer.clear();
            std::getline(stream, stringBuffer);  // Read EOL.
            stream.read(charBuffer, bufferSize);
        }
        return std::shared_ptr<ParticleSpecCatalog>(new ParticleSpecCatalog(specs));
    }
    
    std::ostream& operator << (std::ostream& stream, const ParticleSpecCatalog& catalog)
    {
        using map_specs_t = ParticleSpecCatalog::map_specs_t;
        using iter_t = map_specs_t::const_iterator;
        using pair_t = std::pair<std::string, particle_spec_ptr_t>;
        
        map_specs_t specs = catalog.specs_;
        std::size_t size = specs.size();
        std::size_t index = 0;
        for (iter_t iter = specs.begin(); iter != specs.end(); ++iter ) {
            pair_t pair = *iter;
            auto spec = *pair.second;
            stream << spec;
            if ( index < (size - 1) ) {
                stream << std::endl;
            }
            index += 1;
        }
        return stream;
    }
}
