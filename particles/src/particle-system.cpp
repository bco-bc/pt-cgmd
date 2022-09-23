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
#include "simploce/conf/p-conf.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <boost/algorithm/string/trim.hpp>

namespace simploce {

    void
    ParticleSystem::validate(const simploce::p_system_ptr_t &particleSystem) {
        // For all pairs of particles.
        util::Logger logger("simploce::ParticleSystem::validate()");
        logger.trace("Entering.");
        logger.info("Validating particle system...");
        auto& all = particleSystem->all_;
        auto numberOfParticles = all.size();
        auto numberOfPairs = numberOfParticles * (numberOfParticles - 1) / 2;
        logger.debug(std::to_string(numberOfParticles) + ": Number of particles.");
        logger.debug(std::to_string(numberOfPairs) + ": Number of particle pairs.");
        std::size_t counter = 0;
        for (auto iter_i = all.begin(); iter_i < all.end() - 1; ++iter_i) {
            auto pi = *iter_i;
            auto id_i = pi->id();
            for (auto iter_j = iter_i + 1; iter_j < all.end(); ++iter_j) {
                auto pj = *iter_j;
                auto id_j = pj->id();
                if ( pi == pj) {
                    std::string msg =
                            "Particles (" + util::toString(id_i) + ", " + util::toString(id_j) +
                            ") : Particles point to same particle.";
                    util::logAndThrow(logger, msg);
                }
                if ( id_i == id_j) {
                    std::string msg =
                            "Particles (" + std::to_string(pi->index()) + ", " + std::to_string(pj->index()) +
                            ": Identical particle identifiers.";
                    util::logAndThrow(logger, msg);
                }
                counter += 1;
                if (counter % 5000000 == 0 && counter > 1) {
                    auto percentage = real_t(counter) / real_t(numberOfPairs) * 100.0;
                    logger.debug(std::to_string(percentage) + ": Percentage of particle pairs validated.");
                }
            }
        }
        logger.info("Done. No errors encountered.");
        logger.trace("Leaving.");
    }

    ParticleSystem::ParticleSystem(ParticleSystem&& particleSystem) noexcept:
        all_{}, mapAll_{}, free_{}, displaceables_(), groups_{}, ids_{}, box_{} {
        all_ = std::move(particleSystem.all_);
        mapAll_ = std::move(particleSystem.mapAll_);
        free_ = std::move(particleSystem.free_);
        displaceables_ = std::move(particleSystem.displaceables_);
        groups_ = std::move(particleSystem.groups_);
        ids_ = std::move(particleSystem.ids_);
        box_ = std::move(particleSystem.box_);
    }

    ParticleSystem&
    ParticleSystem::operator = (ParticleSystem&& particleSystem) noexcept {
        all_ = std::move(particleSystem.all_);
        mapAll_ = std::move(particleSystem.mapAll_);
        free_ = std::move(particleSystem.free_);
        displaceables_ = std::move(particleSystem.displaceables_);
        groups_ = std::move(particleSystem.groups_);
        ids_ = std::move(particleSystem.ids_);
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
        ids_.emplace(id);
        return particle;
    }

    p_ptr_t
    ParticleSystem::addParticle(const simploce::id_t &particleId,
                                const std::string &name,
                                const simploce::spec_ptr_t &spec) {
        int index = int(this->numberOfParticles());   // Starts with 0.
        auto particle = this->createParticle_(particleId, index, name, spec);
        if (spec->isFree()) {
            this->addFree(particle);
        } else {
            this->add(particle);
        }
        ids_.emplace(particleId);
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
    ParticleSystem::numberOfDisplaceableParticles() const {
        return displaceables_.size();
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
        auto iter = mapAll_.find(id);
        return iter != mapAll_.end() ? iter->second : nullptr;
        //return util::find(id, all_);
    }

    bool
    ParticleSystem::contains(const id_t& id) const {
        return ids_.find(id) != ids_.end();
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
    ParticleSystem::setOriginToCenterOfMass() {
        auto cm = properties::centerOfMass(all_);
        for (auto& p : all_) {
            auto r = p->position();
            r -= cm;
            p->position(r);
        }
    }

    ParticleSystem::ParticleSystem() : all_{}, mapAll_{}, displaceables_{}, groups_{}, ids_{}, box_{} {
    }

    void
    ParticleSystem::add(const p_ptr_t& particle)
    {
        util::Logger logger{"simploce::ParticleSystem::add()"};
        logger.trace("Entering.");
        //logger.debug(util::toString(particle->id()) + ": Identifier particle being added.");
        /*
        if ( this->contains(particle->id()) ) {
            std::string msg =
                    util::toString(particle->id()) +
                    ": Particle identifier already in use in this particle system.";
            util::logAndThrow(logger, msg);
        }
        auto iter = std::find(all_.begin(), all_.end(), particle);
        if ( iter != all_.end() ) {
            std::string msg =
                    util::toString(particle->id()) +
                    ": Already added to particle system.";
            util::logAndThrow(logger, msg);
        }
         */
        all_.emplace_back(particle);
        auto pair = std::make_pair(particle->id(), particle);
        mapAll_.emplace(pair);
        displaceables_.emplace_back(particle);

        logger.trace("Leaving.");
    }

    void
    ParticleSystem::remove(const p_ptr_t& particle) {
        auto id = particle->id();
        auto iter = std::find(free_.begin(), free_.end(), particle);
        if ( iter != free_.end() ) {
            free_.erase(iter);
        }
        iter = std::find(displaceables_.begin(), displaceables_.end(), particle);
        if ( iter != displaceables_.end() ) {
            displaceables_.erase(iter);
        }
        mapAll_.erase(id);
        iter = std::find(all_.begin(), all_.end(), particle);
        if ( iter != all_.end()  ) {
            all_.erase(iter);
        }
        auto id_iter = std::find(ids_.begin(), ids_.end(), id);
        if ( id_iter != ids_.end()) {
            ids_.erase(id_iter);
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
    ParticleSystem::addFree(const p_ptr_t& particle)
    {
        util::Logger logger{"simploce::ParticleSystem::addFree()"};
        /*
        auto iter = std::find(free_.begin(), free_.end(), particle);
        if ( iter != free_.end() ) {
            util::logAndThrow(logger,
                              util::toString(particle->id()) + ": Already added to particle system.");
        }
        */
        this->add(particle);
        free_.emplace_back(particle);
    }

    void
    ParticleSystem::addParticleGroup(const pg_ptr_t& particleGroup)
    {
        util::Logger logger{"simploce::ParticleSystem::addParticleGroup"};
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
    ParticleSystem::removeFromDisplaceables(const p_ptr_t& particle) {
        auto iter = std::find(displaceables_.begin(), displaceables_.end(), particle);
        if (iter != displaceables_.end() ) {
            displaceables_.erase(iter);
        }
    }

    void
    ParticleSystem::removeFromDisplaceables(const spec_ptr_t& spec) {
        util::Logger logger("simploce::ParticleSystem::removeFromDisplaceables()");
        logger.trace("Entering.");
        std::vector<p_ptr_t> remove{};
        for (const auto& p: displaceables_) {
            if (p->spec()->name() == spec->name()) {
                remove.emplace_back(p);
            }
        }
        logger.debug(std::to_string(remove.size()) + ": Number of particles to be removed from displaceables.");
        std::size_t counter = 0;
        for (const auto& p: remove) {
            this->removeFromDisplaceables(p);
            counter += 1;
            if ( counter % 1000 == 0 && counter > 1) {
                logger.debug(std::to_string(counter) + ": Number of particles removed.");
            }
        }
        logger.debug(std::to_string(counter) + ": Number of particles removed.");
        logger.trace("Leaving");
    }

    void
    ParticleSystem::readFreeAndGroups(std::istream& stream)
    {
        util::Logger logger("simploce::ParticleSystem::readFreeAndGroups()");
        logger.trace("Entering");

        std::string stringBuffer;
        char buffer[1000];

        // Clear free particles and groups to avoid double counting.
        free_.clear();
        groups_.clear();

        // Read free particles.
        std::size_t nFree;
        stream >> nFree;
        logger.debug(std::to_string(nFree) + ": Number of free particles.");
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

            if (counter % 10000 == 0) {
                logger.debug(std::to_string(counter) + ": Number of free particles read.");
            }
        }
        if (nFree > 0 ) {
            std::getline(stream, stringBuffer);  // Read EOL.
        }

        // Read particle groups.
        std::size_t nGroups;
        stream >> nGroups;
        logger.debug(std::to_string(nGroups) + ": Number of particle groups.");
        std::getline(stream, stringBuffer);  // Read EOL.
        for (std::size_t counter = 0; counter != nGroups; ++counter) {

            // Read constituting particles.
            std::vector<p_ptr_t> particles;
            std::size_t nParticles;
            stream >> nParticles;
            logger.debug(std::to_string(nParticles) + ": Number of particles in current particle group.");
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
            logger.debug(std::to_string(nBonds)+ ": Number of bonds in current particle group.");
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

            if ( counter % 10000 == 0) {
                logger.debug(std::to_string(counter) + ": Number of particle group read.");
            }

        }
        logger.trace("Leaving");
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

    std::vector<p_ptr_t>&
    ParticleSystem::displaceables() {
        return displaceables_;
    }

    void
    ParticleSystem::clear() {
        groups_.clear();
        displaceables_.clear();
        free_.clear();
        all_.clear();
        //mapAll_.clear();
        ids_.clear();
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
        util::Logger logger{"simploce::ParticleSystem::parse()"};
        logger.trace("Entering");

        int nParticles, protonatable;
        stream >> nParticles >> protonatable;
        logger.debug(std::to_string(nParticles) + ": Number of particles.");
        if ( protonatable == conf::PROTONATABLE ) {
            logger.warn("Particle system may be protonatable.");
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
            auto particle = this->addParticle(id, name, spec);

            // Position, velocity.
            particle->readState(stream);

            // Done, next.
            std::getline(stream, stringBuffer);  // Read EOL.

            if ( i % 10000 == 0) {
                logger.debug(std::to_string(i) + ": Number of particles read.");
            }
        }

        // Free atoms and atom groups.
        this->readFreeAndGroups(stream);
        box_ptr_t box = factory::box(0.0);
        stream >> *box;
        this->box(box);
        logger.debug("[" + std::to_string(box->lengthX()) + ", " + std::to_string(box->lengthY()) +
                    std::to_string(box->lengthZ()) + ": Box dimensions.");

        // Log some details.
        logger.debug(std::to_string(this->numberOfParticles()) + ": Number of particles.");
        logger.debug(std::to_string(this->numberOfFreeParticles()) + ": Number of free particles.");
        logger.debug(std::to_string(this->numberOfParticleGroups()) + ": Number of particle groups.");

        logger.trace("Leaving.");
    }

    void
    ParticleSystem::resetSpec_(const simploce::p_ptr_t &particle, const simploce::spec_ptr_t &spec) {
        particle->resetSpec(spec);
    }

    std::ostream&
    operator << (std::ostream& stream, const ParticleSystem& particleSystem) {
        particleSystem.write(stream);
        return stream;
    }

}
