/*
 * File:   physical-system.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 12:05 PM
 */

#ifndef PARTICLE_SYSTEM_HPP
#define PARTICLE_SYSTEM_HPP

#include "p-types.hpp"

namespace simploce {

    /**
     * An entity of interest observed in nature. For instance, gases, liquids, 
     * solids and plasmas as well as molecules, atoms, nuclei, and hadrons. 
     * Consists of interacting 'particles'. A particle system consists of
     * free particles (e.g., an ion) and particle groups (protein residue, molecule, etc.).
     * There are also 'displaceable' particles, which refer to particles that can be, e.g.,
     * moved ('displaced') to a new position. Any particle is assumed to be
     * displaceable unless explicitly removed from the list of displaceable particles.
     * @tparam P Particle type.
     * @tparam PG Particle group type.
     */
    class ParticleSystem {
    public:
        
        // Noncopyable.
        ParticleSystem(const ParticleSystem&) = delete;
        ParticleSystem& operator = (const ParticleSystem&) = delete;
        
        // Movable.
        ParticleSystem(ParticleSystem&& particleSystem) noexcept;
        ParticleSystem& operator = (ParticleSystem&& particleSystem) noexcept;
        
        /**
         * Adds a new particle to this particle system. All arguments are required.
         * @param name Name (does not need to be unique).
         * @param spec Particle specification.
         * @return Added particle.
         */
        p_ptr_t addParticle(const std::string& name,
                            const spec_ptr_t& spec);

        /**
         * Adds a particle group with bonds to this physical system. All arguments are
         * required.
         * @param particles Particles forming a particle group. Each of these particles
         * must already be present in this particle system.
         * @param bonds Bonds between particles, given as pairs of particle identifiers.
         * @return Particle group added.
         */
        pg_ptr_t addParticleGroup(const std::vector<p_ptr_t>& particles,
                                  const std::vector<id_pair_t>& bonds);

        /**
         * Returns total number of particles.
         * @return Number.
         */
        std::size_t numberOfParticles() const;
        
        /**
         * Returns number of free particles, that is the number of particles that are -not- part of
         * any particle group.
         * @return Number.
         */
        std::size_t numberOfFreeParticles() const;
        
        /**
         * Returns number of particle groups.
         * @return Number.
         */
        std::size_t numberOfParticleGroups() const;

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
        charge_t charge() const;
        
        /**
         * Returns total mass.
         * @return Value.
         */
        mass_t mass() const;
        
        /**
         * Resets forces on all particles to zero.
         */
        void resetForces();        
        
        /**
         * Are there any particles in this system?
         * @return Result.
         */
        bool empty() const;
        
        /**
         * Finds particle with given identifier.
         * @param id Identifier.
         * @return Particle, or nullptr if not found.
         */
        p_ptr_t find(const id_t& id) const;
        
        /**
         * is there a particle with with the given identifier.
         * @param id Particle identifier.
         * @return Result.
         */
        bool contains(const id_t& id) const;

        /**
         * Returns particle box.
         * @return Box.
         */
        box_ptr_t box() const;
        
        /**
         * Sets the box.
         * @param box
         */
        void box(box_ptr_t box);

        /**
         * Represents this particle model a protonatable system.
         * @return Result. Default is false.
         */
        virtual bool isProtonatable() const;

        /**
         * Performs a "task" with -all- particles. The given task must expose the operator
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
         * Performs a "task" with all "displaceable" particles. Displaceable means here that the particle
         * in question may undergo a state transition (e.g., new position). Note that this method restricts the
         * task to those particles that can be displaced in contrast to the other 'doWith*' methods.
         * The given task must expose the operator
         * <code>
         * R operator () (const std::vector<p_ptr_t>& displaceables);
         * </code>
         * where 'displaceables' represents all displaceable particles.
         * @tparam R Return type. May be 'void'.
         * @tparam TASK Task type.
         * @param task Task. This may be a lambda expression.
         * @return Result.
         */
        template <typename R, typename TASK>
        R doWithDisplaceables(const TASK& task) { return task(displaceables_); }
        
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
         * Writes this particle system to an output stream. Displaceable particles are -not-
         * stored in the output stream.
         * @param stream Output stream.
         */
        void write(std::ostream& stream) const;
        
        /**
         * Writes the current state to an output stream. The state of -all- particles is
         * written to the output stream.
         * @param stream Output stream.
         */
        void writeState(std::ostream& stream) const;

        /**
         * Read the current state from an input stream. The state of -all- particles is
         * read from the input stream.
         * @param stream Input stream.
         */
        void readState(std::istream& stream);
        
    protected:
        
        /**
         * Constructor. No particles, no box.
         */
        ParticleSystem() : all_{}, free_{}, displaceables_{}, groups_{}, box_{} {}
        
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
         * Reassigns a sequential index to all particles.
         */
        void resetIndex();
        
        /**
         * Adds a free particle. Must -not- already be present.
         * @param particle Free particle.
         */
        void addFree(const p_ptr_t& particle);

        /**
         * Adds particle group. Must -not- already be present. Any particle in this group
         * currently in 'free' is removed from 'free'.
         * @param particleGroup Particle group.
         */
        void addParticleGroup(const pg_ptr_t& particleGroup);

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
         * Removes a particle from the list of displaceable particles.
         * @param particle Particle.
         */
        void removeFromDisplaceables(const p_ptr_t& particle);

        /**
         * Removes particles of given particle specification from the list of
         * displaceable particles.
         * @param spec Particle specification.
         */
        void removeFromDisplaceables(const spec_ptr_t& spec);
        
        /**
         * Reads free particles and particle groups from an input stream.
         * @param stream Input stream.
         */
        void readFreeAndGroups(std::istream& stream);

        /**
         * Returns all particle groups.
         * @return Groups. May be empty.
         */
        std::vector<pg_ptr_t> &groups();
        
        /**
         * Returns all particles.
         * @return Particles. Should not be empty.
         */
        std::vector<p_ptr_t> &all();

        /**
         * Returns free particles.
         * @return Particles. May be empty.
         */
        std::vector<p_ptr_t> &free();

        /**
         * Returns all displaceable particles.
         * @return Particles. May be empty.
         */
        std::vector<p_ptr_t> &displaceables();

        /**
         * Clears this particle system. That is, removes all particles and particle groups.
         */
        void clear();

        /**
         * Generate unique particle identifier.
         * @return Identifier.
         */
        id_t generateParticleId() const;

        /**
         * Assigns particle identifier to given particle.
         * @param id Identifier.
         * @param particle Particle.
         */
        static void assignParticleId(const id_t& id, p_ptr_t& particle);

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

        // Displaceable particles. Particles are always displaceables unless removed from this
        // list.
        std::vector<p_ptr_t> displaceables_;
        
        // Particles groups.
        std::vector<pg_ptr_t> groups_;

        // Particle box.
        box_ptr_t box_;
    };

    /**
     * Write particle system to an output stream.
     * @param stream Output stream.
     * @param particleSystem Particle system.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const ParticleSystem& particleSystem);


}

#endif /* PARTICLE_SYSTEM_HPP */

