/*
 * Author: André H. Juffer.
 * Created on 16/11/2021, 09:48.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/p-properties.hpp"
#include "simploce/particle/p-util.hpp"
#include "simploce/particle/p-conf.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <boost/algorithm/string/trim.hpp>

namespace simploce {

    ParticleSystem::ParticleSystem(ParticleSystem&& particleSystem) noexcept:
        all_{}, free_{}, groups_{}, box_{} {
        all_ = std::move(particleSystem.all_);
        free_ = std::move(particleSystem.free_);
        groups_ = std::move(particleSystem.groups_);
        box_ = std::move(particleSystem.box_);
    }

    ParticleSystem&
    ParticleSystem::operator = (ParticleSystem&& particleSystem) noexcept {
        all_ = std::move(particleSystem.all_);
        free_ = std::move(particleSystem.free_);
        groups_ = std::move(particleSystem.groups_);
        box_ = std::move(particleSystem.box_);
        return *this;
    }

    p_ptr_t
    ParticleSystem::addParticle(const std::string& name,
                                      const spec_ptr_t& spec) {
        id_t id = this->generateParticleId();
        int index = int(this->numberOfParticles());   // Starts with 0.
        auto particle = this->createParticle_(id, index, name, spec);
        if ( spec->isFree() ) {
            this->addFree(particle);
        } else {
            this->add(particle);
        }
        return particle;
    }

    pg_ptr_t
    ParticleSystem::addParticleGroup(const std::vector<p_ptr_t>& particles,
                                     const std::vector<id_pair_t>& bonds) {
        auto group = ParticleGroup::make(particles, bonds);
        this->addParticleGroup(group);
        return group;
    }

    std::size_t
    ParticleSystem::numberOfParticles() const {
        return all_.size();
    }

    std::size_t
    ParticleSystem::numberOfFreeParticles() const {
        return free_.size();
    }

    std::size_t
    ParticleSystem::numberOfParticleGroups() const {
        return groups_.size();
    }

    std::size_t
    ParticleSystem::numberOfSpecifications(const spec_ptr_t& spec) const {
        std::size_t counter = 0;
        for (auto& p : all_) {
            if ( p->spec()->name() == spec->name() ) {
                counter += 1;
            }
        }
        return counter;
    }

    charge_t
    ParticleSystem::charge() const {
        return properties::charge<Particle>(all_);
    }

    mass_t
    ParticleSystem::mass() const {
        return properties::mass<Particle>(all_);
    }

    void
    ParticleSystem::resetForces()
    {
        for (auto& p : all_) {
            p->resetForce();
        }
    }

    bool
    ParticleSystem::empty() const {
        return this->numberOfParticles() == 0;
    }

    p_ptr_t ParticleSystem::find(const id_t& id) const {
        return util::find(id, all_);
    }

    bool
    ParticleSystem::contains(const id_t& id) const {
        return this->find(id) != nullptr;
    }

    box_ptr_t
    ParticleSystem::box() const {
        return box_;
    }

    void
    ParticleSystem::box(box_ptr_t box) {
        box_ = std::move(box);
    }

    bool
    ParticleSystem::isProtonatable() const {
        return false;
    }

    void
    ParticleSystem::write(std::ostream& stream) const
    {
        stream << all_.size() << conf::SPACE << this->isProtonatable() << std::endl;
        for (const auto& p: all_) {
            p->write(stream);
            stream << std::endl;
        }
        stream << free_.size() << std::endl;
        for (const auto& p : free_) {
            stream << std::setw(conf::ID_WIDTH) << util::toString(p->id());
        }
        if ( !free_.empty() ) {
            stream << std::endl;
        }
        stream << groups_.size() << std::endl;
        for (auto iter = groups_.begin(); iter != groups_.end(); ++iter) {
            const auto &group = **iter;
            stream << group;
            if ( iter != (groups_.end() - 1) ) {
                stream << std::endl;
            }
        }
        if ( !groups_.empty()) {
            stream << std::endl;
        }
        stream << *box_;
    }

    void
    ParticleSystem::writeState(std::ostream& stream) const
    {
        for (const auto& p: all_) {
            p->writeState(stream);
        }
        stream << std::endl;
    }

    void
    ParticleSystem::readState(std::istream& stream) {
        std::string stringBuffer;
        for (auto& p: all_) {
            p->readState(stream);
        }
        std::getline(stream, stringBuffer);  // Read EOL.
    }

    void
    ParticleSystem::add(const p_ptr_t& particle)
    {
        util::Logger logger{"simploce::ParticleSystem<P,PG>::add()"};

        if ( this->contains(particle->id()) ) {
            std::string msg =
                    util::toString(particle->id()) +
                    ": Particle ID already in use in particle system.";
            util::logAndThrow(logger, msg);
        }
        auto iter = std::find(all_.begin(), all_.end(), particle);
        if ( iter != all_.end() ) {
            std::string msg =
                    util::toString(particle->id()) +
                    ": Already added to particle system.";
            util::logAndThrow(logger, msg);
        }
        all_.push_back(particle);
    }

    void
    ParticleSystem::remove(const p_ptr_t& particle) {
        auto iter = std::find(free_.begin(), free_.end(), particle);
        if ( iter != free_.end() ) {
            free_.erase(iter);
        }
        iter = std::find(all_.begin(), all_.end(), particle);
        if ( iter != all_.end()  ) {
            all_.erase(iter);
        }
    }

    void ParticleSystem::resetIndex() {
        int index = 0;
        for (auto& p : all_) {
            p->index(index);
            index += 1;
        }
    }

    void
    ParticleSystem::addFree(const p_ptr_t& freeParticle)
    {
        util::Logger logger{"simploce::ParticleSystem<P,PG>::addFree()"};
        if ( this->contains(freeParticle->id()) ) {
            util::logAndThrow(logger,
                              util::toString(freeParticle->id()) + ": Already added to particle system.");
        }
        auto iter = std::find(free_.begin(), free_.end(), freeParticle);
        if ( iter != free_.end() ) {
            util::logAndThrow(logger,
                              util::toString(freeParticle->id()) + ": Already added to particle system.");
        }
        this->add(freeParticle);
        free_.emplace_back(freeParticle);
    }

    void
    ParticleSystem::addParticleGroup(const pg_ptr_t& particleGroup)
    {
        util::Logger logger{"simploce::ParticleSystem<P,PG>::addParticleGroup"};
        for (const p_ptr_t& p : particleGroup->particles() ) {
            if ( !this->contains(p->id()) ) {
                util::logAndThrow(
                        logger,
                        util::toString(p->id()) +
                                                ": Particle in particle group not yet added to particle system."
                );
            }
        }
        for (const auto& g : groups_ ) {
            if (g == particleGroup) {
                util::logAndThrow(
                        logger,
                        util::toString(particleGroup->id()) +
                                                ": Particle group already added to particle system.");
            }
        }
        groups_.emplace_back(particleGroup);
        this->removeFromFree(particleGroup);
    }

    void
    ParticleSystem::removeGroup(const pg_ptr_t& particleGroup) {
        auto iter = std::find(groups_.begin(), groups_.end(), particleGroup);
        if ( iter != groups_.end() ) {
            groups_.erase(iter);
        }
    }

    void ParticleSystem::removeFromFree(const pg_ptr_t& particleGroup) {
        for ( auto& p : particleGroup->particles() ) {
            auto iter = std::find(free_.begin(), free_.end(), p);
            if ( iter != free_.end() ) {
                free_.erase(iter);
            }
        }
    }

    void
    ParticleSystem::readFreeAndGroups(std::istream& stream)
    {
        util::Logger logger("simploce::ParticleSystem::readFreeAndGroups");
        std::string stringBuffer;
        char buffer[1000];

        // Clear free particles and groups to avoid double counting.
        free_.clear();
        groups_.clear();

        // Read free particles.
        std::size_t nFree;
        stream >> nFree;
        std::getline(stream, stringBuffer);  // Read EOL.
        for (std::size_t counter = 0; counter != nFree; ++counter) {

            // Identifier
            stream.read(buffer, conf::ID_WIDTH);
            std::string str(buffer, conf::ID_WIDTH);
            boost::trim(str);
            id_t id = util::toId(str);

            // Find particle.
            p_ptr_t particle = this->find(id);
            if ( particle == nullptr ) {
                util::logAndThrow(logger, util::toString(id) + ": No such free particle.");
            }

            free_.push_back(particle);
        }
        if (nFree > 0 ) {
            std::getline(stream, stringBuffer);  // Read EOL.
        }

        // Read particle groups.
        std::size_t nGroups;
        stream >> nGroups;
        std::getline(stream, stringBuffer);  // Read EOL.
        for (std::size_t counter = 0; counter != nGroups; ++counter) {

            // Read constituting particles.
            std::vector<p_ptr_t> particles;
            std::size_t nParticles;
            stream >> nParticles;
            std::getline(stream, stringBuffer);  // Read EOL.
            for (std::size_t j = 0; j != nParticles; ++j) {

                // Identifier
                stream.read(buffer, conf::ID_WIDTH);
                std::string str(buffer, conf::ID_WIDTH);
                boost::trim(str);
                id_t id = util::toId(str);

                // Find particle.
                p_ptr_t particle = this->find(id);
                if ( particle == nullptr ) {
                    util::logAndThrow(logger, util::toString(id) + ": No such group particle.");
                }

                particles.push_back(particle);
            }
            std::getline(stream, stringBuffer);  // Read EOL.

            // Read bonds.
            std::vector<id_pair_t> bonds;
            std::size_t nBonds;
            stream >> nBonds;
            std::getline(stream, stringBuffer);  // Read EOL.
            for (std::size_t j = 0; j != nBonds; ++j) {
                id_t id1, id2;
                stream >> id1 >> id2;
                if ( !this->contains(id1) || !this->contains(id2) ) {
                    util::logAndThrow(
                            logger,
                            util::toString(id1) + ", " + util::toString(id2) + ": No such particle(s)."
                    );
                }
                id_pair_t bond = std::make_pair(id1, id2);
                bonds.push_back(bond);
                std::getline(stream, stringBuffer);  // Read EOL.
            }

            // Create the group.
            auto group = ParticleGroup::make(particles, bonds);
            this->addParticleGroup(group);
        }
    }

    std::vector<pg_ptr_t>&
            ParticleSystem::groups() {
        return groups_;
    }

    std::vector<p_ptr_t>&
    ParticleSystem::all() {
        return all_;
    }

    std::vector<p_ptr_t>&
    ParticleSystem::free() {
        return free_;
    }

    void
    ParticleSystem::clear() {
        groups_.clear();
        free_.clear();
        all_.clear();
        box_ = factory::box(length_t{0.0});
    }

    id_t
    ParticleSystem::generateParticleId() const {
        id_t id = util::generateId();
        while ( this->contains(id) ) {
            id = util::generateId();
        }
        return id;
    }

    void
    ParticleSystem::assignParticleId(const id_t& id, p_ptr_t& particle) {
        particle->id(id);
    }

    void
    ParticleSystem::parse(std::istream &stream, const spec_catalog_ptr_t &catalog) {
        util::Logger logger{"simploce::ParticleSystem<P,PG>::parse()"};

        int nParticles, protonatable;
        stream >> nParticles >> protonatable;
        if ( protonatable == conf::PROTONATABLE ) {
            logger.debug("Particle model may be protonatable.");
        }
        std::string stringBuffer;
        std::getline(stream, stringBuffer);  // Read EOL.

        // Read particles.
        char buffer[1000];
        for (int i = 0; i != nParticles; ++i) {

            // Particle name
            stream.read(buffer, conf::NAME_WIDTH);
            std::string name(buffer, conf::NAME_WIDTH);
            boost::trim(name);

            // Particle sequence index.
            int index;
            stream >> index;

            // Particle specification.
            stream.read(buffer, conf::NAME_WIDTH);
            std::string specName(buffer, conf::NAME_WIDTH);
            boost::trim(specName);
            auto spec = catalog->lookup(specName);

            // Identifier.
            stream.read(buffer, conf::ID_WIDTH);
            std::string str(buffer, conf::ID_WIDTH);
            boost::trim(str);
            id_t id = util::toId(str);

            // Add particle.
            auto particle = this->addParticle(name, spec);
            this->assignParticleId(id, particle);

            // Position, velocity.
            particle->readState(stream);

            // Done, next.
            std::getline(stream, stringBuffer);  // Read EOL.
        }

        // Free atoms and atom groups.
        this->readFreeAndGroups(stream);
        box_ptr_t box = factory::box(0.0);
        stream >> *box;
        this->box(box);

        // Log some details.
        logger.debug("Number of particles: " + util::toString(this->numberOfParticles()));
        logger.debug("Number of free particles: " + util::toString(this->numberOfFreeParticles()));
        logger.debug("Number of particle groups: " + util::toString(this->numberOfParticleGroups()));

    }

    std::ostream&
    operator << (std::ostream& stream, const ParticleSystem& particleSystem) {
        particleSystem.write(stream);
        return stream;
    }

}
