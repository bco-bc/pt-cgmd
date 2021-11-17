/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on November 16, 2021.
 */

#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/bond.hpp"
#include "simploce/particle/p-properties.hpp"
#include "simploce/particle/p-types.hpp"
#include "simploce/util/util.hpp"
#include <boost/algorithm/string.hpp>
#include <memory>
#include <stdexcept>
#include <utility>

namespace simploce {

    ParticleGroup::ParticleGroup(p_ptr_cont_t  particles) :
            id_{++ID_}, particles_{std::move(particles)}, bonds_{}
    {
        this->validate_();
    }

    ParticleGroup::ParticleGroup(p_ptr_cont_t particles,
                                 bond_cont_t  bonds) :
            id_{++ID_}, particles_{std::move(particles)}, bonds_{std::move(bonds)}
    {
        this->validate_();
    }

    ParticleGroup::ParticleGroup(const std::vector<p_ptr_t>& particles,
                                 const std::vector<id_pair_t>& bonds) :
            id_{++ID_}, particles_{}, bonds_{}
    {
        for (const auto& p : particles) {
            particles_.insert(p);
        }
        for (const auto& pair : bonds) {
            auto p1 = this->find_(pair.first);
            auto p2 = this->find_(pair.second);
            if ( p1 == nullptr || p2 == nullptr ) {
                throw std::domain_error("Particle in bond not in particle group.");
            }
            Bond bond = Bond::makeBond(p1, p2);
            bonds_.push_back(bond);
        }
        this->validate_();
    }

    bool
    ParticleGroup::operator == (const ParticleGroup& group) const {
        return id_ == group.id_;
    }

    std::size_t
    ParticleGroup::id() const {
        return id_;
    }

    bool
    ParticleGroup::contains(const p_ptr_t& p) const {
        return particles_.find(p) != particles_.end();
    }

    charge_t
    ParticleGroup::charge() const {
        return properties::charge(particles_);
    }

    mass_t
    ParticleGroup::mass() const {
        return properties::mass(particles_);
    }

    position_t
    ParticleGroup::position() const {
        return properties::centerOfMass(particles_);
    }

    const ParticleGroup:: p_ptr_cont_t&
    ParticleGroup::particles() const {
        return particles_;
    }

    const ParticleGroup::bond_cont_t&
    ParticleGroup::bonds() const {
        return bonds_;
    }

    pg_ptr_t
    ParticleGroup::make(const std::vector<p_ptr_t>& particles,
                        const std::vector<id_pair_t>& bonds)
    {
        return std::make_shared<ParticleGroup>(particles, bonds);
    }

    p_ptr_t
    ParticleGroup::find_(const id_t& id) const {
        return properties::find(id, particles_);
    }

    void
    ParticleGroup::validate_() const
    {
        if (particles_.empty() ) {
            throw std::domain_error(
                    "ParticleGroup: Must consist of at least two particles."
            );
        }
        for (const auto& p : particles_) {
            if ( p == nullptr ) {
                throw std::domain_error(
                        "ParticleGroup: Missing particle (\"nullptr\")."
                );
            }
        }
        for (const auto& bond : bonds_) {
            if ( !this->contains(bond.getParticleOne()) ||
                 !this->contains(bond.getParticleTwo()) ) {
                throw std::domain_error(
                        "Particle in bond of particle group is not in the containing particle group."
                );
            }
        }
    }

    std::size_t
    ParticleGroup::ID_ = 0;

    std::ostream&
    operator << (std::ostream& stream,
                 const ParticleGroup& group)
    {
        stream << group.particles().size() << std::endl;
        for (const auto& p : group.particles() ) {
            stream << std::setw(conf::ID_WIDTH) << util::toString(p->id());
        }
        stream << std::endl;
        stream << group.bonds().size() << std::endl;
        for (const auto& bond : group.bonds() ) {
            std::string id1 = util::toString(bond.getParticleOne()->id());
            std::string id2 = util::toString(bond.getParticleTwo()->id());
            stream << std::setw(conf::ID_WIDTH) << id1 << std::setw(conf::ID_WIDTH) << id2;
        }

        return stream;
    }

}

