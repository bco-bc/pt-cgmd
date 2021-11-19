/*
 * File:   particle.cpp
 * Author: André H. Juffer, Biocenter Oulu
 *
 * Created on July 30, 2019, 4:06 PM
 */

#include "simploce/particle/p-conf.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/util.hpp"
#include <stdexcept>
#include <iomanip>
#include <utility>

namespace simploce {
    
    Particle::Particle(id_t id,
                       std::size_t index,
                       std::string name,
                       spec_ptr_t spec) :
        id_{std::move(id)}, index_{index}, name_{std::move(name)}, spec_{std::move(spec)},
        r_{}, v_{}, f_{} {
        if ( id_.empty() ) {
            throw std::domain_error("A particle identifier must be provided.");
        }
        if ( name_.empty() ) {
            throw std::domain_error("A particle name must be provided.");
        }
        if (!spec_) {
            throw std::domain_error("A particle specification must be provided.");
        }
    }
    
    Particle::~Particle() = default;

    bool
    Particle::operator == (const Particle& p) const {
        return id_ == p.id_;
    }
    
    id_t
    Particle::id() const {
        return id_;
    }
    
    std::size_t 
    Particle::index() const {
        return index_;
    }
    
    std::string 
    Particle::name() const {
        return name_;
    }
    
    spec_ptr_t 
    Particle::spec() const {
        return spec_;
    }
    
    charge_t 
    Particle::charge() const {
        return spec_->charge();
    }
    
    mass_t 
    Particle::mass() const {
        return spec_->mass();
    }
    
    position_t
    Particle::position() const {
        return r_;
    }
    
    void 
    Particle::position(const position_t& r) {
        r_ = r; 
    }
    
    velocity_t
    Particle::velocity() const {
        return v_;
    }
    
    void 
    Particle::velocity(const velocity_t& v) {
        v_ = v;
    }
    
    force_t
    Particle::force() const {
        return f_; 
    }
    
    void 
    Particle::force(const force_t& f) {
        f_ = f;
    }
    
    void 
    Particle::resetForce() {
        f_ = force_t{0.0, 0.0, 0.0};
    }

    bool
    Particle::isIon() const {
        return spec_->isIon();
    }
    
    void 
    Particle::write(std::ostream& stream) const {
        const auto space = conf::SPACE;

        stream << std::setw(conf::NAME_WIDTH) << this->name();
        stream << std::setw(conf::INTEGER_WIDTH) << this->index();
        stream << std::setw(conf::NAME_WIDTH) << this->spec()->name();
        stream << std::setw(conf::ID_WIDTH) << util::toString(this->id());
        stream << space << this->position();
        stream << space << this->velocity();
    }
    
    void 
    Particle::writeState(std::ostream& stream) const {
        stream << this->position() << this->velocity();
    }
    
    void 
    Particle::readState(std::istream& stream) {
        real_t x, y, z, vx, vy, vz;
        stream >> x >> y >> z >> vx >> vy >> vz;
        position_t r{x, y, z};
        velocity_t v{vx, vy, vz};
        this->position(r);
        this->velocity(v);
    }

    void
    Particle::id(const id_t& id) {
        id_ = id;
    }
    
    void 
    Particle::resetSpec(const spec_ptr_t &spec) {
        if ( !spec ) {
            throw std::domain_error("A particle specification must be provided.");
        }
        spec_ = spec;
    }
    
    std::ostream& 
    operator << (std::ostream& stream, 
                 const Particle& particle) {
        particle.write(stream);
        return stream;
    }
    
}
