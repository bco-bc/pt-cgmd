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
 * File:   physical-system.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 12:05 PM
 */

#ifndef PARTICLE_MODEL_HPP
#define PARTICLE_MODEL_HPP

#include "particle-group.hpp"
#include "pproperties.hpp"
#include "ptypes.hpp"
#include "pconf.hpp"
#include "simploce/util/util.hpp"
#include <memory>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <utility>

namespace simploce {
    
    /**
     * An entity of interest observed in nature. For instance, gases, liquids, 
     * solids and plasmas as well as molecules, atoms, nuclei, and hadrons. 
     * Consists of interacting 'particles'. This class is noncopyable, but it is 
     * movable.
     * @param P Particle type.
     */
    template <typename P, typename PG>
    class ParticleModel {
    public:
        
        // Noncopyable.
        ParticleModel(const ParticleModel&) = delete;
        ParticleModel& operator = (const ParticleModel&) = delete;
        
        // Movable.
        ParticleModel(ParticleModel&& particleModel);
        ParticleModel& operator = (ParticleModel&& particleModel);
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Particle group pointer type.
         */
        using pg_ptr_t = std::shared_ptr<PG>;
        
        /**
         * Returns number of particles.
         * @return Number.
         */
        std::size_t numberOfParticles() const { return all_.size(); }
        
        /**
         * Returns of free particles.
         * @return Number.
         */
        std::size_t numberOfFreeParticles() const { return free_.size(); }
        
        /**
         * Returns number of particle groups.
         * @return Number.
         */
        std::size_t numberOfParticleGroups() const { return groups_.size(); }
        
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
         * Resets forces on particles to zero.
         */
        void resetForces();        
        
        /**
         * Returns number of protonation sites.
         * @return Number.
         */
        virtual std::size_t numberOfProtonationSites() const = 0;
        
        /**
         * Returns total number of bound protons.
         * @return Number, always >= 0.
         */
        virtual std::size_t protonationState() const = 0;
        
        /**
         * Are there any particles in this physical system?
         * @return Result.
         */
        bool empty() const { return this->numberOfParticles() == 0; }
        
        /**
         * Finds particle with given identifier.
         * @param id Identifier.
         * @return Particle, or nullptr if not found.
         */
        p_ptr_t find(std::size_t id) const { return properties::find(id, all_); }
        
        /**
         * Contains a particle with given identifier.
         * @param id Particle identifier.
         * @return Result.
         */
        bool contains(std::size_t id) const { return this->find(id) != nullptr; }
        
        /**
         * Performs a "task" with all particles. The given task must expose the 
         * operator 
         * <code>
         * R operator () (const std::vector<p_ptr_t>& all);
         * </code>
         * where 'all' represents all the particles in this model.
         * @param task Task of type TASK. This may be a lambda expression.
         * @return Result of type R. May be void.
         */
        template <typename R, typename TASK>
        R doWithAll(const TASK& task) { return task(all_); }
        
        /**
         * Performs a "task" with all and free particles, as well as the 
         * particle groups. The given task must expose the operator
         * <code>
         *  R operator () (const std::vector<p_ptr_t>& all,
         *                 const std::vector<p_ptr_t>& free,
         *                 const std::vector<pg_ptr_t>& groups);
         * </code>
         * where 'all' represents -all- particles, 'free' refers to those particles
         * not in any particle group (e.g. ions), and 'groups' are all particle 
         * groups, respectively.
         * @param task Task of type TASK. This may be a lambda expression.
         * @return Result of type R. May be void.
         */
        template <typename R, typename TASK>
        R doWithAllFreeGroups(const TASK& task) { return task(all_, free_, groups_); }
        
        /**
         * Writes this particle model to an output stream.
         * @param stream Output stream.
         */
        void write(std::ostream& stream) const;
        
        /**
         * Writes the current state to an output stream.
         * @param stream Output stream.
         */
        virtual void writeState(std::ostream& stream) const;
        
    protected:
        
        /**
         * Constructor. No particles.
         */
        ParticleModel() : all_{}, free_{}, groups_{} {}
        
        /**
         * Adds particle.
         * @param p Particle.
         */
        void add(const p_ptr_t& p);
        
        /**
         * Adds free particle.
         * @param particle Free particle.
         */
        void addFree(const p_ptr_t& fp);
        
        /**
         * Adds particle group.
         * @param g Particle group.
         */
        void addGroup(const pg_ptr_t& pg);
        
        /**
         * Reads free particles and particle groups from input stream.
         * @param stream Input stream.
         */
        void readFreeAndGroups(std::istream& stream);
        
        /**
         * Returns groups.
         * @return 
         */
        std::vector<pg_ptr_t>& groups() { return groups_; }
        
        /**
         * Return particles.
         */
        std::vector<p_ptr_t>& all() { return all_; }
                
    private:
        
        // All particles
        std::vector<p_ptr_t> all_;
        
        // Free particles, not in any group.
        std::vector<p_ptr_t> free_;
        
        // Particles groups.
        std::vector<pg_ptr_t> groups_;
    };
    
    template <typename P, typename PG>
    ParticleModel<P,PG>::ParticleModel(ParticleModel&& pm) :
        all_{}, free_{}, groups_{}
    {
        all_ = std::move(pm.all_);
        free_ = std::move(pm.free_);
        groups_ = std::move(pm.groups_);
    }
        
    template <typename P, typename PG>
    ParticleModel<P,PG>& ParticleModel<P,PG>::operator = (ParticleModel&& pm)
    {
        all_ = std::move(pm.all_);
        free_ = std::move(pm.free_);
        groups_ = std::move(pm.groups_);
        return *this;
    }
    
    template <typename P, typename PG>
    void 
    ParticleModel<P,PG>::resetForces()
    {
        std::for_each(all_.begin(), all_.end(), [] (P& p) { 
            p->resetForce(); 
        });
    }
    
    template <typename P, typename PG>
    void 
    ParticleModel<P,PG>::write(std::ostream& stream) const
    {
        const char space = conf::SPACE;
        
        stream << all_.size() << std::endl;
        for (const auto& p: all_) {
            p->write(stream);
            stream << std::endl;            
        }
        stream << free_.size() << std::endl;        
        for (const auto& p : free_) {
            stream << space << p->id();
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
    }
    
    template <typename P, typename PG>
    void 
    ParticleModel<P,PG>::writeState(std::ostream& stream) const
    {
        for (const auto& p: all_) {
            p->writeState(stream);        
        }
    }
    
    template <typename P, typename PG>
    void 
    ParticleModel<P,PG>::add(const p_ptr_t& p)
    {
        if ( this->contains(p->id()) ) {
            std::string msg = 
                util::toString(p->id()) +  ": Already added to particle model."; 
            throw std::domain_error(msg);
        }
        all_.push_back(p);
    }
    
    template <typename P, typename PG>
    void 
    ParticleModel<P,PG>::addFree(const p_ptr_t& fp) 
    {
        if ( this->contains(fp->id()) ) {
            throw std::domain_error(
                util::toString(fp->id()) +  ": Already added to particle model."   
            );
        }
        this->add(fp);
        free_.push_back(fp); 
    }
    
    template <typename P, typename PG>
    void
    ParticleModel<P,PG>::addGroup(const pg_ptr_t& pg)
    {
        for (const p_ptr_t& p : pg->particles() ) {
            if ( !this->contains(p->id()) ) {
                throw std::domain_error(
                    util::toString(p->id()) + ": Not yet added to particle model." 
                );
            }
        }
        for (auto g : groups_ ) {
            if ( g == pg) {
                throw std::domain_error(
                    "Particle group already added to particle model."
                );
            }
        }
        groups_.push_back(pg);
    }
        
    template <typename P, typename PG>
    void
    ParticleModel<P,PG>::readFreeAndGroups(std::istream& stream)
    {
        std::string stringBuffer;
        
        // Read free particles.
        std::size_t nfree;
        stream >> nfree;
        std::getline(stream, stringBuffer);  // Read EOL.
        for (std::size_t counter = 0; counter != nfree; ++counter) {
            std::size_t id;
            stream >> id;
            p_ptr_t particle = this->find(id);
            if ( particle == nullptr ) {
                throw std::domain_error(util::toString(id) + ": No such free particle.");
            }
            free_.push_back(particle);
        }
        if ( nfree > 0 ) {
            std::getline(stream, stringBuffer);  // Read EOL.
        }
        
        // Read particle groups.
        std::size_t ngroups;
        stream >> ngroups;
        std::getline(stream, stringBuffer);  // Read EOL.
        for ( std::size_t counter = 0; counter != ngroups; ++counter) {
            
            // Read constituting particles.
            std::vector<p_ptr_t> particles;
            std::size_t nparticles;
            stream >> nparticles;
            std::getline(stream, stringBuffer);  // Read EOL.
            for (std::size_t j = 0; j != nparticles; ++j) {
                std::size_t id;
                stream >> id;
                p_ptr_t particle = this->find(id);
                if ( particle == nullptr ) {
                    throw std::domain_error(util::toString(id) + ": No such particle.");
                }
                particles.push_back(particle);
            }                        
            std::getline(stream, stringBuffer);  // Read EOL.
            
            // Read bonds.
            std::vector<id_pair_t> bonds;
            std::size_t nbonds;
            stream >> nbonds;
            std::getline(stream, stringBuffer);  // Read EOL.
            for (std::size_t j = 0; j != nbonds; ++j) {
                std::size_t id1, id2;
                stream >> id1 >> id2;
                id_pair_t bond = std::make_pair(id1, id2);
                bonds.push_back(bond);
                std::getline(stream, stringBuffer);  // Read EOL.
            }
            
            // Create the group.
            pg_ptr_t group = PG::make(particles, bonds);
            this->addGroup(group);            
        }        
    }
      
}

#endif /* PARTICLE_MODEL_HPP */

