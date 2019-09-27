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
    
    spec_ptr_t 
    ParticleSpecCatalog::lookup(const std::string& name) const
    {
        auto iter = specs_.find(name);
        if ( iter == specs_.end() ) {
            throw std::domain_error(
                name + ": No specification available for this particle."
            );
        }
        return (*iter).second;
    }
    
    spec_ptr_t 
    ParticleSpecCatalog::molecularWater() const
    {
        return this->lookup("mH2O");
    }
    
    spec_catalog_ptr_t 
    ParticleSpecCatalog::create(std::istream& stream)
    {
        map_specs_t specs{};
        
        const size_t bufferSize = conf::NAME_WIDTH;
        char charBuffer[bufferSize];
        std::string stringBuffer;
        
        // Skip header (first line).
        std::getline(stream, stringBuffer);
        
        bool protonatable;
        stream >> protonatable;  // Represents a protonatable particle?
        while ( stream.good() ) {
            stream.read(charBuffer, bufferSize);  // Extra spaces.
            stream.read(charBuffer, bufferSize);  // Should hold spec name.
            std::string name = std::string(charBuffer, bufferSize); 
            boost::trim(name);
            real_t charge, mass, radius;        
            stream >> mass >> charge >> radius;            
            if ( !protonatable ) {
                spec_ptr_t spec = 
                    ParticleSpec::create(name, charge, mass, radius);
                auto pair = std::make_pair(name, spec);
                specs.insert(pair);
            } else {
                real_t pKa;
                bool continuous;
                stream >> pKa >> continuous;
                spec_ptr_t spec = 
                    ParticleSpec::create(name, charge, mass, radius, pKa, continuous);
                auto p = std::make_pair(name, spec);
                specs.insert(p);                
            }
            
            stringBuffer.clear();
            std::getline(stream, stringBuffer);  // Read EOL.
            
            // Next.
            stream >> protonatable;
        }
        return spec_catalog_ptr_t(new ParticleSpecCatalog{specs});
    }
    
    void 
    ParticleSpecCatalog::write(std::ostream& stream) const
    {
        std::size_t counter = 0;
        for (auto iter = specs_.begin(); iter != specs_.end(); ++iter)
        {
            counter += 1;
            auto pair = *iter;
            auto spec = *pair.second;
            stream << spec;
            if ( counter < specs_.size() ) {
                stream << std::endl;
            }
        }
    }
    
    std::ostream& 
    operator << (std::ostream& stream, 
                 const ParticleSpecCatalog& catalog)
    {
        catalog.write(stream);
        return stream;
    }
}
