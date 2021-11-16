/*
 * File:   particle-spec.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 2:01 PM
 */

#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/p-conf.hpp"
#include <stdexcept>
#include <utility>

namespace simploce {

    spec_ptr_t
    ParticleSpec::create(const std::string& name,
                         const charge_t& charge,
                         const mass_t& mass,
                         const radius_t& radius,
                         bool free,
                         const std::string& description)
    {
        return ParticleSpec::create_(name, charge, mass, radius, 0.0, false, free, description);
    }

    spec_ptr_t
    ParticleSpec::create(const std::string& name,
                         const charge_t& charge,
                         const mass_t& mass,
                         const radius_t& radius,
                         const pKa_t& pKa,
                         bool free,
                         const std::string& description)
    {
        return ParticleSpec::create_(name, charge, mass, radius, pKa, true, free, description);
    }

    spec_ptr_t
    ParticleSpec::createFrom(const spec_ptr_t& spec,
                             const std::string& name,
                             charge_t charge,
                             bool free,
                             const std::string& description)
    {
        return ParticleSpec::create_(name,
                                     charge,
                                     spec->mass(),
                                     spec->radius(),
                                     spec->pKa(),
                                     spec->isProtonatable(),
                                     free,
                                     description);
    }

    bool
    ParticleSpec::operator == (const ParticleSpec& particleSpec) {
        return name_ == particleSpec.name();
    }

    std::string
    ParticleSpec::name() const
    {
        return name_;
    }
    
    charge_t 
    ParticleSpec::charge() const
    {
        return charge_;
    }
    
    mass_t 
    ParticleSpec::mass() const
    {
        return mass_;
    }
    
    radius_t 
    ParticleSpec::radius() const
    {
        return radius_;
    }
    
    pKa_t 
    ParticleSpec::pKa() const
    {
        return pKa_;
    }
    
    bool 
    ParticleSpec::isProtonatable() const
    {
        return protonatable_;
    }
    
    bool
    ParticleSpec::isIon() {
        if ( name_ == "Na+" ||
             name_ == "Cl-") {
            return true;
        }
        return false;
    }

    bool
    ParticleSpec::isFree() {
        return free_;
    }

    std::string
    ParticleSpec::description() const {
        return description_;
    }
    
    void
    ParticleSpec::write(std::ostream& stream) const
    {
        stream.setf(std::ios::scientific);
        stream.precision(conf::PRECISION);
        stream << std::setw(conf::NAME_WIDTH) << name_;
        stream << conf::SPACE << protonatable_;
        stream << conf::SPACE << free_;
        stream << std::setw(conf::REAL_WIDTH) << mass_;
        stream << std::setw(conf::REAL_WIDTH) << charge_;
        stream << std::setw(conf::REAL_WIDTH) << radius_;
        stream << std::setw(conf::REAL_WIDTH) << pKa_;
        stream << conf::SPACE << description_;
    }

    spec_ptr_t
    ParticleSpec::create_(const std::string& name,
                          const charge_t& charge,
                          const mass_t& mass,
                          const radius_t& radius,
                          const pKa_t& pKa,
                          bool protonatable,
                          bool free,
                          const std::string& description) {
        if ( name.empty() ) {
            throw std::domain_error("An unique particle specification name must be provided.");
        }
        if ( mass() < 0.0) {
            throw std::domain_error("A mass must be a non-negative number.");
        }
        if ( radius() < 0.0 ) {
            throw std::domain_error("A radius must be a non-negative number.");
        }
        if ( description.empty() ) {
            throw std::domain_error("A particle specification description must be provided.");
        }
        return spec_ptr_t(
            new ParticleSpec(name, charge, mass, radius, pKa, protonatable, free, description)
        );
    }

    
    ParticleSpec::ParticleSpec(std::string name,
                               charge_t charge,
                               mass_t mass,
                               radius_t radius,
                               pKa_t pKa, 
                               bool protonatable,
                               bool free,
                               std::string description) :
        name_{std::move(name)}, charge_{charge}, mass_{mass}, radius_{radius}, pKa_{pKa},
        protonatable_{protonatable}, free_{free}, description_(std::move(description)) {
    }
        
    std::ostream& 
    operator << (std::ostream& stream, 
                 const ParticleSpec& spec)
    {
        spec.write(stream);
        return stream;
    }    
}
