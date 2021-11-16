/*
 * File:   particle.hpp
 * Author: André H. Juffer, Biocenter Oulu
 *
 * Created on July 30, 2019, 4:05 PM
 */

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "p-types.hpp"
#include <memory>
#include <string>

namespace simploce {
    
    /**
     * A recognizable and identifiable unit composing a physical system at any time. Has location,
     * feels forces, and may be moving. Also, has charge and mass.
     */
    class Particle {
    public:
        
        // Noncopyable.
        Particle(const Particle&) = delete;
        Particle& operator = (const Particle&) = delete;

        // Compares identifiers.
        bool operator == (const Particle& p) const;
        
        virtual ~Particle();
        
        /**
         * Returns unique particle identifier. This uniquely identifies a given 
         * particle in a given particle collection. Once assigned, the particle's ID will never
         * change.
         * @return Identifier.
         */
        id_t id() const;
        
        /**
         * Returns particle sequential index. Must -not- be used as an particle identifier.
         * @return Index, always >= 0. This index may change, e.g., if particles are removed.
         * @see #id()
         */
        std::size_t index() const;
        
        /**
         * Returns particle name. May not be unique.
         * @return Name.
         */
        std::string name() const;
        
        /**
         * Returns particle specification. Multiple particles may have the same specification.
         * @return Specification.
         */
        spec_ptr_t spec() const;
        
        /**
         * Returns particle charge.
         * @return Charge.
         */
        virtual charge_t charge() const;
        
        /**
         * Returns particle mass.
         * @return Mass. Its value is always > 0.
         */
        virtual mass_t mass() const;
        
        /**
         * Returns current position.
         * @return Position.
         */
        position_t position() const;
        
        /**
         * Sets position.
         * @param r New position.
         */
        void position(const position_t& r);
        
        /**
         * Returns current velocity.
         * @return Velocity.
         */
        velocity_t velocity() const;
        
        /**
         * Sets velocity.
         * @param v Velocity.
         */
        void velocity(const velocity_t& v);
        
        /**
         * Returns force acting on this particle.
         * @return Force.
         */
        force_t force() const;
        
        /**
         * Sets force acting on this particle.
         * @param f Force.
         */
        void force(const force_t& f);
        
        /**
         * Resets force acting on this particle to zero.
         */
        void resetForce();

        /**
         * Is this particle an ion?
         * @return Result.
         */
        bool isIon() const;
        
        /**
         * Writes this particle to an output stream.
         * @param stream Output stream.
         */
        virtual void write(std::ostream& stream) const;
        
        /**
         * Writes particle state (position, velocity) to an output stream.
         * @param stream Output stream.
         */
        virtual void writeState(std::ostream& stream) const;
        
        /**
         * Reads particle state (position, velocity) from an input stream.
         * @param stream Input stream.
         */
        virtual void readState(std::istream& stream);
        
    protected:
        
        /**
         * Constructor. All arguments must be provided.
         * @param id Unique particle identifier.
         * @param index Unique particle index.
         * @param name Particle name. Does not need to be unique.
         * @param spec Particle specification.
         */
        Particle(id_t id,
                 std::size_t index, 
                 std::string name,
                 spec_ptr_t spec);

        /**
         * Assigns particle sequential index.
         * @param index Index.
         */
        void index(size_t index) { index_ = index; }
                        
    private:
        
        template <typename P, typename PG>
        friend class ParticleSystem;

        template <typename P>
        friend class ProtonationSite;

        // Sets particle id.
        void id(const id_t& id);
        
        // Assigns particle specification.
        void resetSpec(const spec_ptr_t& spec);

        id_t id_;
        std::size_t index_;
        std::string name_;
        spec_ptr_t spec_;
        position_t r_;
        velocity_t v_;
        force_t f_;
        
    };
    
    /**
     * Writes particle to output stream.
     * @param stream Output stream.
     * @param particle Particle.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const Particle& particle);
}



#endif /* PARTICLE_HPP */

