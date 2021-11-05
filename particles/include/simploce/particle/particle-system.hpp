/*
 * File:   physical-system.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 12:05 PM
 */

#ifndef PARTICLE_SYSTEM_HPP
#define PARTICLE_SYSTEM_HPP

#include "particle-group.hpp"
#include "particle-spec.hpp"
#include "particle-spec-catalog.hpp"
#include "p-properties.hpp"
#include "p-types.hpp"
#include "p-conf.hpp"
#include "p-factory.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <boost/algorithm/string.hpp>
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <utility>

namespace simploce {

    /**
     * An entity of interest observed in nature. For instance, gases, liquids, 
     * solids and plasmas as well as molecules, atoms, nuclei, and hadrons. 
     * Consists of interacting 'particles'. A particle system consists of
     * free particles (e.g., an ion) and particle groups (protein residue, molecule, etc.).
     * @tparam P Particle type.
     * @tparam PG Particle group type.
     */
    template <typename P, typename PG>
    class ParticleSystem {
    public:
        
        // Noncopyable.
        ParticleSystem(const ParticleSystem&) = delete;
        ParticleSystem& operator = (const ParticleSystem&) = delete;
        
        // Movable.
        ParticleSystem(ParticleSystem&& particleSystem) noexcept;
        ParticleSystem& operator = (ParticleSystem&& particleSystem) noexcept;
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Particle group pointer type.
         */
        using pg_ptr_t = std::shared_ptr<PG>;

        /**
         * Addsa new particle to this particle system. All arguments are required.
         * @param name Name (does not need to be unique).
         * @param spec Particle specification.
         * @return Added particle.
         */
        p_ptr_t addParticle(const std::string& name,
                            const spec_ptr_t& spec);

        /**
         * Returns total number of particles.
         * @return Number.
         */
        std::size_t numberOfParticles() const { return all_.size(); }
        
        /**
         * Returns number of free particles, that is the number of particles that are -not- part of
         * any particle group.
         * @return Number.
         */
        std::size_t numberOfFreeParticles() const { return free_.size(); }
        
        /**
         * Returns number of particle groups.
         * @return Number.
         */
        std::size_t numberOfParticleGroups() const { return groups_.size(); }

        /**
         * Returns number of particles with given specification.
         * @param spec Particle specification.
         * @return Number.
         */
        std::size_t numberOfSpecifications(const spec_ptr_t& spec) const;
        
        /**
         * Returns total charge.
         * @return Value.
         */
        charge_t charge() const { return properties::charge<P>(all_); }
        
        /**
         * Returns total mass.
         * @return Value.
         */
        mass_t mass() const { return properties::mass<P>(all_); }
        
        /**
         * Resets forces on all particles to zero.
         */
        void resetForces();        
        
        /**
         * Are there any particles in this system?
         * @return Result.
         */
        bool empty() const { return this->numberOfParticles() == 0; }
        
        /**
         * Finds particle with given identifier.
         * @param id Identifier.
         * @return Particle, or nullptr if not found.
         */
        p_ptr_t find(const id_t& id) const { return properties::find(id, all_); }
        
        /**
         * is there a particle with with the given identifier.
         * @param id Particle identifier.
         * @return Result.
         */
        bool contains(const id_t& id) const { return this->find(id) != nullptr; }

        /**
         * Returns particle box.
         * @return Box.
         */
        box_ptr_t box() const { return box_; }
        
        /**
         * Sets the box.
         * @param box
         */
        void box(const box_ptr_t& box) { box_ = factory::box(box->size()); }

        /**
         * Represents this particle model a protonatable system.
         * @return Result. Default is false.
         */
        virtual bool isProtonatable() const { return false; }

        /**
         * Performs a "task" with -all- particles. The given task must expose the
         * operator 
         * <code>
         * R operator () (const std::vector<p_ptr_t>& all);
         * </code>
         * where 'all' represents all the particles in this model.
         * @tparam R Return type. May be 'void'.
         * @tparam TASK Task type.
         * @param task Task. This may be a lambda expression.
         * @return Result.
         */
        template <typename R, typename TASK>
        R doWithAll(const TASK& task) { return task(all_); }
        
        /**
         * Performs a "task" with all, free particles, as well as with the
         * particle groups. The given task must expose the operator
         * <code>
         *  R operator () (std::vector<p_ptr_t>& all,
         *                 std::vector<p_ptr_t>& free,
         *                 std::vector<pg_ptr_t>& groups);
         * </code>
         * where 'all' represents -all- particles, 'free' refers to those particles
         * not in any particle group (e.g. ions), and 'groups' are all particle 
         * groups, respectively.
         * @tparam R Return type. May be 'void'.
         * @tparam TASK Task type.
         * @param task Task. This may be a lambda expression.
         * @return Result.
         */
        template <typename R, typename TASK>
        R doWithAllFreeGroups(const TASK& task) { return task(all_, free_, groups_); }

        /**
         * Writes this particle system to an output stream.
         * @param stream Output stream.
         */
        void write(std::ostream& stream) const;
        
        /**
         * Writes the current state to an output stream.
         * @param stream Output stream.
         */
        void writeState(std::ostream& stream) const;
        
    protected:
        
        /**
         * Constructor. No particles.
         */
        ParticleSystem() : all_{}, free_{}, groups_{}, box_{factory::box(length_t{0})} {}
        
        /**
         * Adds particle.
         * @param particle Particle.
         */
        void add(const p_ptr_t& particle);

        /**
         * Completely removes particle from this particle system.
         * @param particle Particle to be removed.
         */
        void remove(const p_ptr_t& particle);

        /**
         * Reassigns sequential index to all particles.
         */
        void resetIndex();
        
        /**
         * Adds a free particle. Must -not- already be present.
         * @param freeParticle Free particle.
         */
        void addFree(const p_ptr_t& freeParticle);
        
        /**
         * Adds particle group. Must -not- already be present. Any particle in this group
         * currently in 'free' is removed from 'free'.
         * @param particleGroup Particle group.
         */
        void addGroup(const pg_ptr_t& particleGroup);

        /**
         * Removes a particle group.
         * @param particleGroup Particle group.
         */
        void removeGroup(const pg_ptr_t& particleGroup);

        /**
         * Removes particle in group from free particles.
         * @param particleGroup Particle group.
         */
        void removeFromFree(const pg_ptr_t& particleGroup);
        
        /**
         * Reads free particles and particle groups from an input stream.
         * @param stream Input stream.
         */
        void readFreeAndGroups(std::istream& stream);

        /**
         * Returns all particle groups.
         */
        std::vector<pg_ptr_t>& groups() { return groups_; }
        
        /**
         * Returns all particles.
         */
        std::vector<p_ptr_t>& all() { return all_; }

        /**
         * Returns free particles.
         */
        std::vector<p_ptr_t>& free() { return free_; }

        /**
         * Clears this particle system. That is, removes all particles and particle groups.
         */
        void clear();

        /**
         * Generate unique particle identifier.
         * @return Identifier.
         */
        id_t generateParticleId();

        /**
         * Assigns particle identifier to given particle.
         * @param id Identifier.
         * @param particle Particle.
         */
        void assignParticleId(const id_t& id, p_ptr_t& particle);

        /**
         * Obtains particle system from an input stream.
         * @param stream Input stream.
         * @param catalog Particle specification catalog.
         */
        void parse(std::istream& stream,
                   const spec_catalog_ptr_t& catalog);

    private:

        /**
         * Creates particle.
         * @param id Unique particle identifier.
         * @param index Particle sequence number.
         * @param name Name (does not need to be unique).
         * @param r Position.
         * @param spec Particle specification.
         * @return Created particle.
         */
        virtual p_ptr_t createParticle_(const id_t& id,
                                        int index,
                                        const std::string& name,
                                        const spec_ptr_t& spec) = 0;

        // All particles
        std::vector<p_ptr_t> all_;
        
        // Free particles, not in any group.
        std::vector<p_ptr_t> free_;
        
        // Particles groups.
        std::vector<pg_ptr_t> groups_;

        // Particle box.
        box_ptr_t box_;
    };

    template <typename P, typename PG>
    ParticleSystem<P,PG>::ParticleSystem(ParticleSystem&& particleSystem) noexcept:
        all_{}, free_{}, groups_{}, box_{} {
        all_ = std::move(particleSystem.all_);
        free_ = std::move(particleSystem.free_);
        groups_ = std::move(particleSystem.groups_);
        box_ = std::move(particleSystem.box_);
    }
        
    template <typename P, typename PG>
    ParticleSystem<P,PG>& ParticleSystem<P,PG>::operator = (ParticleSystem&& particleSystem) noexcept {
        all_ = std::move(particleSystem.all_);
        free_ = std::move(particleSystem.free_);
        groups_ = std::move(particleSystem.groups_);
        box_ = std::move(particleSystem.box_);
        return *this;
    }

    template <typename P, typename PG>
    typename ParticleSystem<P,PG>::p_ptr_t
    ParticleSystem<P,PG>::addParticle(const std::string& name,
                                      const spec_ptr_t& spec) {
        id_t id = this->generateParticleId();
        int index = int(this->numberOfParticles()) + 1;
        auto particle = this->createParticle_(id, index, name, spec);
        if ( spec->isIon() ) {
            this->addFree(particle);
        } else {
            this->add(particle);
        }
        return particle;
    }


    template <typename P, typename PG>
    std::size_t ParticleSystem<P,PG>::numberOfSpecifications(const spec_ptr_t& spec) const {
        std::size_t counter = 0;
        for (auto& p : all_) {
            if ( p->spec()->name() == spec->name() ) {
                counter += 1;
            }
        }
        return counter;
    }
    
    template <typename P, typename PG>
    void 
    ParticleSystem<P,PG>::resetForces()
    {
        for (auto& p : all_) {
            p->resetForce();
        }
    }
    
    template <typename P, typename PG>
    void 
    ParticleSystem<P,PG>::write(std::ostream& stream) const
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
            const PG& group = **iter;
            stream << group;
            if ( iter != (groups_.end() - 1) ) {
                stream << std::endl;
            }
        }
        if ( groups_.size() != 0) {
            stream << std::endl;
        }
        stream << *box_;
    }
    
    template <typename P, typename PG>
    void 
    ParticleSystem<P,PG>::writeState(std::ostream& stream) const
    {
        for (const auto& p: all_) {
            p->writeState(stream);
        }
    }

    template <typename P, typename PG>
    void 
    ParticleSystem<P,PG>::add(const p_ptr_t& particle)
    {
        util::Logger logger{"simploce::ParticleSystem<P,PG>::add()"};

        if ( this->contains(particle->id()) ) {
            std::string msg =
                    util::toString(particle->id()) + ": Particle ID already in use in particle system.";
            util::logAndThrow(logger, msg);
        }
        auto iter = std::find(all_.begin(), all_.end(), particle);
        if ( iter != all_.end() ) {
            std::string msg =
                    util::toString(particle->id()) + ": Already added to particle system.";
            util::logAndThrow(logger, msg);
        }
        all_.push_back(particle);
    }

    template <typename P, typename PG>
    void
    ParticleSystem<P,PG>::remove(const p_ptr_t& particle) {
        using iter_t = typename std::vector<p_ptr_t>::iterator;

        iter_t iter = std::find(free_.begin(), free_.end(), particle);
        if ( iter != free_.end() ) {
            free_.erase(iter);
        }
        iter = std::find(all_.begin(), all_.end(), particle);
        if ( iter != all_.end()  ) {
            all_.erase(iter);
        }
    }

    template <typename P, typename PG>
    void ParticleSystem<P,PG>::resetIndex() {
        int index = 1;
        for (auto& p : all_) {
            p->index(index);
            index += 1;
        }
    }
    
    template <typename P, typename PG>
    void 
    ParticleSystem<P,PG>::addFree(const p_ptr_t& freeParticle)
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
        free_.push_back(freeParticle);
    }
    
    template <typename P, typename PG>
    void
    ParticleSystem<P,PG>::addGroup(const pg_ptr_t& particleGroup)
    {
        util::Logger logger{"simploce::ParticleSystem<P,PG>::addGroup"};
        for (const p_ptr_t& p : particleGroup->particles() ) {
            if ( !this->contains(p->id()) ) {
                util::logAndThrow(
                        logger,
                        util::toString(p->id()) + ": Particle in particle group not yet added to particle system."
                );
           }
        }
        for (auto g : groups_ ) {
            if (g == particleGroup) {
                util::logAndThrow(
                        logger,
                        util::toString(particleGroup->id()) + ": Particle group already added to particle system.");
            }
        }
        groups_.push_back(particleGroup);
        this->removeFromFree(particleGroup);
    }

    template <typename P, typename PG>
    void
    ParticleSystem<P,PG>::removeGroup(const pg_ptr_t& particleGroup) {
        using iter_t = typename std::vector<pg_ptr_t>::iterator;
        iter_t iter = std::find(groups_.begin(), groups_.end(), particleGroup);
        if ( iter != groups_.end() ) {
            groups_.erase(iter);
        }
    }

    template <typename P, typename PG>
    void ParticleSystem<P,PG>::removeFromFree(const pg_ptr_t& particleGroup) {
        for ( auto& p : particleGroup->particles() ) {
            auto iter = std::find(free_.begin(), free_.end(), p);
            if ( iter != free_.end() ) {
                free_.erase(iter);
            }
        }
    }
        
    template <typename P, typename PG>
    void
    ParticleSystem<P,PG>::readFreeAndGroups(std::istream& stream)
    {
        util::Logger logger("simploce::ParticleSystem<P,PG>::readFreeAndGroups");
        std::string stringBuffer;
        char buffer[1000];

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
            pg_ptr_t group = PG::make(particles, bonds);
            this->addGroup(group);            
        }
    }

    template <typename P, typename PG>
    void
    ParticleSystem<P,PG>::clear() {
        groups_.clear();
        free_.clear();
        all_.clear();
        box_ = factory::box(length_t{0.0});
    }

    template <typename P, typename PG>
    id_t
    ParticleSystem<P,PG>::generateParticleId() {
        static std::set<id_t> ids;
        id_t id = util::generateId();
        while ( ids.find(id) != ids.end() ) {
            id = util::generateId();
        }
        ids.insert(id);
        return id;
    }

    template <typename P, typename PG>
    void
    ParticleSystem<P,PG>::assignParticleId(const id_t& id, p_ptr_t& particle) {
        particle->id(id);
    }

    template <typename P, typename PG>
    void
    ParticleSystem<P,PG>::parse(std::istream &stream, const spec_catalog_ptr_t &catalog) {
        util::Logger logger{"simploce::ParticleSystem<P,PG>::parse()"};

        int nParticles, protonatable;
        stream >> nParticles >> protonatable;
        if ( protonatable == conf::PROTONATABLE ) {
            logger.debug("Particle model may contain titrating sites.");
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
            real_t x, y, z, vx, vy, vz;

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
      
}

#endif /* PARTICLE_SYSTEM_HPP */

