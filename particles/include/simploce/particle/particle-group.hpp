/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
#include "pproperties.hpp"
#include "ptypes.hpp"
#include "pconf.hpp"
#include <memory>
#include <vector>
#include <set>
#include <stdexcept>
#include <iostream>

namespace simploce {
    
  /**
   * Holds particles that form a logical unit such as a molecule. May include bonds.
   * Interactions between particles in the same group may be handled separately from 
   * interactions between particles in different groups.
   * @param P Particle type.
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
         * Particle container type.
         */
        using p_ptr_cont_t = std::set<p_ptr_t>;
        
        /**
         * Bond container type.
         */
        using bond_cont_t = std::vector<Bond<P>>;
        
        /**
         * Constructor. Particles, no bonds.
         * @param particles Constituting particles.
         */
        ParticleGroup(const p_ptr_cont_t& particles);
        
        /**
         * Constructor. Particles and bonds.
         * @param particles Constituting particles.
         * @param bonds Bonds. Particles forming bonds must be constituting particles.
         */
        ParticleGroup(const p_ptr_cont_t& particles, const bond_cont_t& bonds);
        
        /**
         * Constructor.
         * @param particles Constituting particles.
         * @param pairs Holds identifiers of particles forming bonds. The 
         * particles must be constituting particles.
         */
        ParticleGroup(const std::vector<p_ptr_t>& particles, 
                      const std::vector<id_pair_t>& pairs);
        
        // Noncopyable
        ParticleGroup(const ParticleGroup&) = delete;
        ParticleGroup& operator = (const ParticleGroup&) = delete;
        
        /**
         * Is the given particle in this group.
         * @param particle Particle.
         * @return Result.
         */
        bool contains(const p_ptr_t& p) const { return particles_.find(p) != particles_.end(); }
        
        /**
         * Returns group's charge.
         * @return Charge.
         */
        charge_t charge() const { return properties::charge<P>(particles_); }
        
        /**
         * Returns group's mass.
         * @return Mass.
         */
        mass_t mass() const { return properties::mass<P>(particles_); }
        
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
         * Create particle group.
         * @param particles Constituting particles.
         * @param bonds Holds identifiers of particles forming bonds. The 
         * particles must be constituting particles.
         * @return Particle group.
         */
        static pg_ptr_t make(const std::vector<p_ptr_t>& particles, 
                             const std::vector<id_pair_t>& bonds);
        
    private:
        
        p_ptr_t find_(std::size_t id) const { return properties::find<P>(id, particles_); }
        
        void validate_() const;
        
        p_ptr_cont_t particles_;
        bond_cont_t bonds_;
    };
    
    template <typename P>
    ParticleGroup<P>::ParticleGroup(const p_ptr_cont_t& particles) :
        particles_{particles}, bonds_{}
    {
        this->validate_();
    }
        
    template <typename P>
    ParticleGroup<P>::ParticleGroup(const p_ptr_cont_t& particles,
                                    const bond_cont_t& bonds) :
        particles_{particles}, bonds_{bonds}
    {
        this->validate_();
    }
    
    template <typename P>
    ParticleGroup<P>::ParticleGroup(const std::vector<p_ptr_t>& particles,
                                    const std::vector<id_pair_t>& pairs) :
        particles_{}, bonds_{}
    {
        for (auto p : particles) {
            particles_.insert(p);
        }
        for (auto pair : pairs) {
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
    typename ParticleGroup<P>::pg_ptr_t 
    ParticleGroup<P>::make(const std::vector<p_ptr_t>& particles, 
                           const std::vector<id_pair_t>& bonds)
    {
        return std::make_shared<ParticleGroup<P>>(particles, bonds);
    }    
        
    template <typename P>
    void ParticleGroup<P>::validate_() const
    {
        if (particles_.empty() ) {
            throw std::domain_error(
                "ParticleGroup: Must consist of at least two particles"
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
                    "Particle in bond of particle group is not in the "
                    "containing particle group."
                );
            }
        }
    }
    
    template <typename P>
    std::ostream& operator << (std::ostream& stream, const ParticleGroup<P>& group)
    {        
        const char space = conf::SPACE;
        
        stream << group.particles().size() << std::endl;
        for (auto p : group.particles() ) {
            stream << space << p->id();
        }
        stream << std::endl;
        stream << group.bonds().size() << std::endl;
        for (auto bond : group.bonds() ) {
            stream << space << bond.getParticleOne()->id()
                   << space << bond.getParticleTwo()->id();
        }
        
        return stream;
    }
        
}

#endif /* PARTICLE_GROUP_HPP */

