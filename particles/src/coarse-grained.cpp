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
 * File:   coarse-grained.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:38 PM
 */

#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/protonatable-bead.hpp"
#include "simploce/particle/protonating-bead.hpp"
#include "simploce/particle/pconf.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-group.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <utility>
#include <memory>

namespace simploce {
    
    static std::size_t indexGenerator_(const CoarseGrained& cg)
    {
        return cg.numberOfParticles();
    }
        
    CoarseGrained::CoarseGrained() :
        ParticleModel<Bead,bead_group_t>{}, discrete_{}, continouos_{}
    {        
    }
    
    bead_ptr_t CoarseGrained::addBead(std::size_t id,
                                      const std::string& name, 
                                      const position_t& r, 
                                      const bead_spec_ptr_t& spec)
    {
        auto index = indexGenerator_(*this);
        bead_ptr_t bead = Bead::create(id, index, name, spec);
        bead->position(r);
        this->add(bead);
        return bead;
    }
            
    dprot_bead_ptr_t CoarseGrained::addProtonatableBead(std::size_t id,
                                                       const std::string& name, 
                                                       const position_t& r,
                                                       std::size_t protonationState,
                                                       const bead_spec_ptr_t& spec)
    {  
        auto index = indexGenerator_(*this);
        dprot_bead_ptr_t bead = 
            ProtonatableBead::create(id, index, name, protonationState, spec);
        bead->position(r);
        discrete_.push_back(bead);
        this->add(bead);
        return bead;
    }
    
    cprot_bead_ptr_t CoarseGrained::addProtonatingBead(std::size_t id,
                                                       const std::string& name, 
                                                       const position_t& r,
                                                       std::size_t protonationState,
                                                       const bead_spec_ptr_t& spec)
    {
        auto index = indexGenerator_(*this);
        cprot_bead_ptr_t bead = 
            ProtonatingBead::create(id, index, name, protonationState, spec);
        bead->position(r);
        continouos_.push_back(bead);
        this->add(bead);
        return bead;
    }
    
    bead_group_ptr_t CoarseGrained::addBeadGroup(const std::vector<bead_ptr_t>& beads, 
                                                 const std::vector<id_pair_t>& pairs)
    {
        bead_group_ptr_t group = std::make_shared<ParticleGroup<Bead>>(beads, pairs);
        this->addGroup(group);
        return group;
    }
    
    std::size_t CoarseGrained::numberOfProtonationSites() const
    {
        return discrete_.size() + continouos_.size();
    }
                
    std::size_t CoarseGrained::protonationState() const
    {
        std::size_t numberOfProtons{0};
        for (auto d: discrete_) {
            numberOfProtons += d->protonationState();
        }
        for (auto c: continouos_) {
            numberOfProtons += c->protonationState();
        }
        return numberOfProtons;
    }
    
    cg_ptr_t CoarseGrained::readFrom(std::istream& stream, 
                                     const spec_catalog_ptr_t& catalog)
    {
        cg_ptr_t cg = std::make_shared<CoarseGrained>();
        
        std::string stringBuffer;
        
        std::size_t nbeads;
        stream >> nbeads;
        std::getline(stream, stringBuffer);  // Read EOL.
                
        const size_t bufferSize = conf::NAME_WIDTH;
        char charBuffer[bufferSize];
        
        // Read beads.
        for (std::size_t counter = 0; counter != nbeads; ++counter) {
            stream.read(charBuffer, bufferSize);
            std::string name = std::string(charBuffer, bufferSize);
            boost::trim(name);
            stream.read(charBuffer, bufferSize);
            std::string specName = std::string(charBuffer, bufferSize);
            boost::trim(specName);
            bead_spec_ptr_t spec = catalog->lookup(specName);
            std::size_t id;
            stream >> id;
            real_t x, y, z, vx, vy, vz;
            stream >> x >> y >> z >> vx >> vy >> vz;
            position_t r{x, y, z};
            velocity_t v{vx, vy, vz};
            std::size_t protonatable;
            stream >> protonatable;
            if ( protonatable == conf::PROTONATABLE) {
                // Discrete changes in protonation state.
                std::size_t protonationState;            
                stream >> protonationState;
                dprot_bead_ptr_t bead = 
                    cg->addProtonatableBead(id, name, r, protonationState, spec);
                bead->velocity(v);
            } else if ( protonatable == conf::PROTONATING) {
                // Continuous change in protonation state.
                std::size_t protonationState;
                real_t x, I;
                stream >> protonationState >> x >> I;
                cprot_bead_ptr_t bead =
                    cg->addProtonatingBead(id, name, r, protonationState, spec);
                bead->velocity(v);
                bead->state(x);
                bead->current(I);
            } else {
                // Regular bead.
                bead_ptr_t bead = cg->addBead(id, name, r, spec);
                bead->velocity(v);
            }
            
            std::getline(stream, stringBuffer);  // Read EOL.
        }
                
        cg->readFreeAndGroups(stream);
        
        return cg; 
    }
    
    std::ostream& operator << (std::ostream& stream, const CoarseGrained& cg)
    {
        cg.write(stream);
        return stream;
    }
    
}

