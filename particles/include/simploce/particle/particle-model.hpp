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
#include <memory>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>

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
        bool contains(std::size_t id) const;
        
        /**
         * Performs a "task" with all particles. The given task must expose the 
         * operator 
         * <code>
         * R operator () (const std::vector<p_ptr_t>& all);
         * </code>
         * where 'all' represents all the particles in this model.
         * @param task Task of type TASK. This may also be a lamdba expression.
         * @return Result of type R.
         */
        template <typename R, typename TASK>
        R doWithAll(const TASK task) { return task(all_); }
        
        /**
         * Performs a 'task' with all and free particles, as well as the 
         * particle groups. The given task must expose the operator
         * <code>
         *  R operator () (const std::vector<p_ptr_t>& all,
         *                 const std::vector<p_ptr_t>& free,
         *                 const std::vector<pg_ptr_t>& groups);
         * </code>
         * where 'all' represents -all- particles, 'free' refers to those particles
         * not in any particle group (e.g. ions), and 'groups' are all particle groups.
         * @param task Task of type TASK. This may also be a lamdba expression.
         * @return Result of type R.
         */
        template <typename R, typename TASK>
        R doWithAllFreeGroups(const TASK task) { return task(all_, free_, groups_); }
        
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
    void ParticleModel<P,PG>::resetForces()
    {
        for (p_ptr_t& p : all_) {
            p->resetForce();
        }
    }
    
    template <typename P, typename PG>
    bool ParticleModel<P,PG>::contains(std::size_t id) const
    {
        return this->find(id) != nullptr;
    }
    
    template <typename P, typename PG>
    void ParticleModel<P,PG>::write(std::ostream& stream) const
    {
        const char space = conf::SPACE;
        
        stream << all_.size() << std::endl;
        for (const p_ptr_t& p: all_) {
            p->write(stream);
            stream << std::endl;            
        }
        stream << free_.size() << std::endl;        
        for (const p_ptr_t& p : free_) {
            stream << space << p->index();
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
    void ParticleModel<P,PG>::writeState(std::ostream& stream) const
    {
        const char space = conf::SPACE;
        
        for (const p_ptr_t& p: all_) {
            p->writeState(stream);        
        }
    }
    
    template <typename P, typename PG>
    void ParticleModel<P,PG>::add(const p_ptr_t& p)
    {
        if ( this->contains(p->index()) ) {
            throw std::domain_error("Particle is already in particle model.");
        }
        all_.push_back(p);
    }
    
    template <typename P, typename PG>
    void ParticleModel<P,PG>::addFree(const p_ptr_t& fp) 
    {
        if ( this->contains(fp->index()) ) {
            throw std::domain_error("Particle is already in particle model.");
        }
        all_.push_back(fp);
        free_.push_back(fp); 
    }
    
    template <typename P, typename PG>
    void ParticleModel<P,PG>::addGroup(const pg_ptr_t& pg)
    {
        for (const p_ptr_t& p : pg->particles() ) {
            if ( !this->contains(p->index()) ) {
                throw std::domain_error(
                    "Particle in particle group is not in physical system."
                );
            }
        }
        for (auto g : groups_ ) {
            if ( g == pg) {
                throw std::domain_error("Particle group is already in physical system.");
            }
        }
        groups_.push_back(pg);
    }
        
}

#endif /* PARTICLE_MODEL_HPP */

