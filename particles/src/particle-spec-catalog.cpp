/*
 * File:   particle-catalog.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 8, 2019, 2:29 PM
 */

#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/conf/p-conf.hpp"
#include "simploce/util/logger.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <sstream>
#include <utility>

namespace simploce {

    ParticleSpecCatalog::ParticleSpecCatalog(map_specs_t specs) :
       specs_{std::move(specs)} {
    }

    bool
    ParticleSpecCatalog::hasSpecification(const std::string &name) const {
        try {
            this->lookup(name);
            return true;
        } catch (std::exception& exception) {
            return false;
        }
    }
    
    spec_ptr_t 
    ParticleSpecCatalog::lookup(const std::string& name) const
    {
        util::Logger logger("simploce::ParticleSpecCatalog::lookup(const std::string& name)");
        auto iter = specs_.find(name);
        if ( iter == specs_.end() ) {
            std::string message = name + ": No specification available for this particle.";
            util::logAndThrow(logger, message);
        }
        return (*iter).second;
    }

    spec_ptr_t
    ParticleSpecCatalog::lookupByElementName(const std::string &name) const {
        util::Logger logger("simploce::ParticleSpecCatalog::lookupByElementName(const std::string &name)");
        if (name.empty()) {
            util::logAndThrow(logger, "Name must be provided.");
        }
        // Assume first character specifies element name.
        const char ch = name[0];
        switch (ch) {
            case 'C': return this->C(); break;
            case 'H': return this->H(); break;
            case 'O': return this->O(); break;
            case 'N': return this->N(); break;
            default: return this->lookup(std::string{ch});
        }
    }
    
    spec_ptr_t 
    ParticleSpecCatalog::molecularWater() const
    {
        return this->lookup("mH2O");
    }
    
    spec_ptr_t 
    ParticleSpecCatalog::H() const
    {
        return this->lookup("H");
    }
    
    spec_ptr_t 
    ParticleSpecCatalog::C() const
    {
        return this->lookup("C");
    }
    
    spec_ptr_t 
    ParticleSpecCatalog::O() const
    {
        return this->lookup("O");
    }
    
    spec_ptr_t 
    ParticleSpecCatalog::N() const {
        return this->lookup("N");
    }

    spec_ptr_t
    ParticleSpecCatalog::staticBP() const {
        return this->lookup("SBP");
    }

    spec_catalog_ptr_t
    ParticleSpecCatalog::obtainFrom(std::istream& stream)
    {
        util::Logger logger("simploce::ParticleSpecCatalog::parseIt()");
        map_specs_t specs{};

        char charBuffer[100];
        std::string stringBuffer;
        
        // Skip header (first line).
        std::getline(stream, stringBuffer);

        // Read first specification name.
        stream.read(charBuffer, conf::NAME_WIDTH);  // Should hold spec name.
        while ( stream.good() ) {
            std::string name = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(name);
            real_t charge, mass, radius, pKa;
            bool protonatable, free;
            stream >> protonatable >> free >> mass >> charge >> radius >> pKa;
            stringBuffer.clear();
            std::string description;
            std::getline(stream, description);
            boost::trim(description);
            if ( !protonatable ) {
                auto spec = ParticleSpec::create(name, charge, mass, radius, free, description);
                auto pair = std::make_pair(name, spec);
                if (specs.count(name) != 0) {
                    util::logAndThrow(logger, name + ": Double entry for this specification in catalog.");
                }
                specs.insert(pair);
            } else {
                auto spec = ParticleSpec::create(name, charge, mass, radius, pKa, free, description);
                auto p = std::make_pair(name, spec);
                specs.insert(p);                
            }

            // Read next specification name, if any.
            stringBuffer.clear();
            stream.read(charBuffer, conf::NAME_WIDTH);  // Should hold unique specification name.
        }

        // Log some information.
        std::string nSpecs = boost::lexical_cast<std::string>(specs.size());
        logger.debug("Number of particle specifications: " + nSpecs);

        return spec_catalog_ptr_t(new ParticleSpecCatalog{specs});
    }
    
    void 
    ParticleSpecCatalog::write(std::ostream& stream) const
    {
        std::size_t counter = 0;
        std::string header =
            "      Name  protonatable-bead?  free-particle? Mass-deprotonated  Charge-deprotonated  Radius  pKa  Short Description";
        stream << header << std::endl;
        for (auto iter = specs_.begin(); iter != specs_.end(); ++iter) {
            counter += 1;
            auto pair = *iter;
            auto spec = pair.second;
            stream << *spec;
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
