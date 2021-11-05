/*
 * File:   particle-group.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 15, 2019, 4:53 PM
 */

#ifndef PARTICLE_GROUP_HPP
#define PARTICLE_GROUP_HPP

#include "bond.hpp"
#include "particle.hpp"
#include "particle-group.hpp"
#include "p-properties.hpp"
#include "p-types.hpp"
#include "p-conf.hpp"
#include "simploce/util/util.hpp"
#include <memory>
#include <vector>
#include <set>
#include <stdexcept>
#include <iostream>

namespace simploce {
    
  /**
   * Holds particles that form a logical unit such as a molecule. May include bonds
   * between particles.
   * @tparam P Particle type.
   */
    template <typename P>
    class ParticleGroup {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;

        /**
         * Particle group pointer type.
         */
        using pg_ptr_t = std::shared_ptr<ParticleGroup<P>>;
        
        /**
         * Particle pointer container type.
         */
        using p_ptr_cont_t = std::set<p_ptr_t>;
        
        /**
         * Bonds container type.
         */
        using bond_cont_t = std::vector<Bond<P>>;
        
        /**
         * Constructor. Particles, no bonds.
         * @param particles Constituting particles.
         */
        explicit ParticleGroup(const p_ptr_cont_t& particles);
        
        /**
         * Constructor. Particles and bonds.
         * @param particles Constituting particles.
         * @param bonds Bonds. Particles forming bonds must be constituting particles.
         */
        explicit ParticleGroup(const p_ptr_cont_t& particles,
                               const bond_cont_t& bonds);
        
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
        bool operator == (const ParticleGroup& group) { return id_ == group.id_; }

        /**
         * Returns group identifier.
         * @return Identifier.
         */
        std::size_t id() const { return id_; }
        
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
        virtual charge_t charge() const { return properties::charge<P>(particles_); }
        
        /**
         * Returns group's total mass.
         * @return Mass.
         */
        virtual mass_t mass() const { return properties::mass<P>(particles_); }
        
        /**
         * Returns this group's position.
         */
        position_t position() const { return properties::centerOfMass<P>(particles_); }
        
        /**
         * Returns particles.
         * @return 
         */
        const p_ptr_cont_t& particles() const { return particles_; }
        
        /**
         * Returns bonds.
         * @return Bonds
         */
        const bond_cont_t& bonds() const { return bonds_; }
        
        /**
         * Creates a particle group.
         * @param particles Constituting particles.
         * @param bonds Holds identifiers of particles forming bonds. The 
         * particles must be constituting particles.
         * @return Particle group.
         */
        static pg_ptr_t make(const std::vector<p_ptr_t>& particles, 
                             const std::vector<id_pair_t>& bonds);
        
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
        
        p_ptr_cont_t particles_;
        bond_cont_t bonds_;
    };
    
    template <typename P>
    ParticleGroup<P>::ParticleGroup(const p_ptr_cont_t& particles) :
        id_{++ID_}, particles_{particles}, bonds_{}
    {
        this->validate_();
    }
        
    template <typename P>
    ParticleGroup<P>::ParticleGroup(const p_ptr_cont_t& particles,
                                    const bond_cont_t& bonds) :
        id_{++ID_}, particles_{particles}, bonds_{bonds}
    {
        this->validate_();
    }
    
    template <typename P>
    ParticleGroup<P>::ParticleGroup(const std::vector<p_ptr_t>& particles,
                                    const std::vector<id_pair_t>& bonds) :
        id_{++ID_}, particles_{}, bonds_{}
    {
        for (auto p : particles) {
            particles_.insert(p);
        }
        for (auto pair : bonds) {
            auto p1 = this->find_(pair.first);
            auto p2 = this->find_(pair.second);
            if ( p1 == nullptr || p2 == nullptr ) {
                throw std::domain_error("Particle in bond not in particle group.");
            }
            Bond<P> bond = Bond<P>::makeBond(p1, p2);
            bonds_.push_back(bond);
        }
        this->validate_();
    }
        
    template <typename P> 
    bool 
    ParticleGroup<P>::contains(const p_ptr_t& p) const { 
        return particles_.find(p) != particles_.end();
    }    
    
    template <typename P>    
    typename ParticleGroup<P>::pg_ptr_t 
    ParticleGroup<P>::make(const std::vector<p_ptr_t>& particles, 
                           const std::vector<id_pair_t>& bonds)
    {
        return std::make_shared<ParticleGroup<P>>(particles, bonds);
    }    
        
    template <typename P>
    typename ParticleGroup<P>::p_ptr_t 
    ParticleGroup<P>::find_(const id_t& id) const {
        return properties::find<P>(id, particles_); 
    }
    
    template <typename P>
    void 
    ParticleGroup<P>::validate_() const
    {
        if (particles_.empty() ) {
            throw std::domain_error(
                "ParticleGroup: Must consist of at least two particles."
            );
        }
        for (auto p : particles_) {
            if ( p == nullptr ) {
                throw std::domain_error(
                    "ParticleGroup: Missing particle (\"nullptr\")."
                );
            }
        }
        for (auto bond : bonds_) {
            if ( !this->contains(bond.getParticleOne()) ||
                 !this->contains(bond.getParticleTwo()) ) {
                throw std::domain_error(
                    "Particle in bond of particle group is not in the containing particle group."
                );
            }
        }
    }
    
    template <typename P>
    std::size_t
    ParticleGroup<P>::ID_ = 0;
    
    template <typename P>
    std::ostream& 
    operator << (std::ostream& stream, 
                 const ParticleGroup<P>& group)
    {        
        stream << group.particles().size() << std::endl;
        for (auto p : group.particles() ) {
            stream << std::setw(conf::ID_WIDTH) << util::toString(p->id());
        }
        stream << std::endl;
        stream << group.bonds().size() << std::endl;
        for (auto bond : group.bonds() ) {
            std::string id1 = util::toString(bond.getParticleOne()->id());
            std::string id2 = util::toString(bond.getParticleTwo()->id());
            stream << std::setw(conf::ID_WIDTH) << id1 << std::setw(conf::ID_WIDTH) << id2;
        }
        
        return stream;
    }
        
}

#endif /* PARTICLE_GROUP_HPP */

