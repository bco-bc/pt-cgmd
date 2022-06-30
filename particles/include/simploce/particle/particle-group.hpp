/*
 * File:   particle-group.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 15, 2019, 4:53 PM
 */

#ifndef PARTICLE_GROUP_HPP
#define PARTICLE_GROUP_HPP

#include "p-types.hpp"
#include <set>
#include <vector>
#include <iostream>

namespace simploce {
    
  /**
   * Holds particles that form a logical unit such as a molecule. May specify bonds
   * between two particles, define bond angles for three particles, etc.
   */
    class ParticleGroup {
    public:
        
        /**
         * Particle pointer container type.
         */
        using p_ptr_cont_t = std::vector<p_ptr_t>;
        
        /**
         * Bonds container type.
         */
        using bond_cont_t = std::vector<Bond>;
        
        /**
         * Constructor. Particles, no bonds.
         * @param particles Constituting particles.
         */
        explicit ParticleGroup(p_ptr_cont_t particles);
        
        /**
         * Constructor. Particles and bonds.
         * @param particles Constituting particles.
         * @param bonds Bonds. Particles forming bonds must be constituting particles.
         */
        explicit ParticleGroup(p_ptr_cont_t particles,
                               bond_cont_t bonds);
        
        /**
         * Constructor.
         * @param particles Constituting particles.
         * @param bonds Holds identifiers of particles forming bonds. The
         * particles must be constituting particles.
         */
        ParticleGroup(const std::vector<p_ptr_t>& particles, 
                      const std::vector<id_pair_t>& bonds);
        
        // Noncopyable
        ParticleGroup(const ParticleGroup&) = delete;
        ParticleGroup& operator = (const ParticleGroup&) = delete;
        
        /**
         * Equality operator.
         * @param group Another group.
         * @return Result.
         */
        bool operator == (const ParticleGroup& group) const;

        /**
         * Returns group identifier.
         * @return Identifier.
         */
        std::size_t id() const;
        
        /**
         * Is the given particle in this group.
         * @param particle Particle.
         * @return Result.
         */
        bool contains(const p_ptr_t& p) const;

        /**
         * Returns group's total charge.
         * @return Charge.
         */
        virtual charge_t charge() const;
        
        /**
         * Returns group's total mass.
         * @return Mass.
         */
        virtual mass_t mass() const;
        
        /**
         * Returns this group's position.
         */
        position_t position() const;
        
        /**
         * Returns particles.
         * @return 
         */
        const p_ptr_cont_t& particles() const;
        
        /**
         * Returns bonds.
         * @return Bonds
         */
        const bond_cont_t& bonds() const;
        
        /**
         * Creates a particle group.
         * @param particles Constituting particles.
         * @param bonds Holds identifiers of particles forming bonds. The 
         * particles must be constituting particles.
         * @return Particle group.
         */
        static pg_ptr_t make(const std::vector<p_ptr_t>& particles, 
                             const std::vector<id_pair_t>& bonds);

        /**
         * Returns pairs of particles in this group not forming a bond with each other,
         * are both involved in the same angle definition, etc.
         * @return Particle pairs.
         */
        const std::vector<std::pair<p_ptr_t, p_ptr_t>>& nonBondedParticlePairs();
        
    private:

        static std::size_t ID_;

        std::size_t id_;

        /**
         * Returns particle with given identifier.
         * @param id Particle identifier.
         * @return Particle.
         */
        p_ptr_t find_(const id_t& id) const;
        
        void validate_() const;

        /**
         * Defines pairs of particles in this group not forming a bond with each other,
         * are both involved in the same angle definition, etc.
         */
        void defineNonBonded_();
        
        p_ptr_cont_t particles_;
        bond_cont_t bonds_;
        std::vector<std::pair<p_ptr_t, p_ptr_t>> nonBonded_;
    };
    
    std::ostream& operator << (std::ostream &stream, const ParticleGroup &group);

}

#endif /* PARTICLE_GROUP_HPP */

