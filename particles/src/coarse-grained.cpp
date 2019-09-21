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
#include "simploce/particle/protonation-site.hpp"
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
        return cg.numberOfParticles() + 1;
    }
    
    
    CoarseGrained::CoarseGrained() :
        ParticleModel<Bead,bead_group_t>{}, protonatableBeads_{}
    {        
    }
    
    bead_ptr_t CoarseGrained::addBead(const std::string& name, 
                                      const position_t& r, 
                                      const bead_spec_ptr_t& spec)
    {
        auto index = indexGenerator_(*this);
        bead_ptr_t bead = Bead::create(index, name, spec);
        bead->position(r);
        this->add(bead);
        return bead;
    }
            
    prot_bead_ptr_t CoarseGrained::addProtonatableBead(const std::string& name, 
                                                       const position_t& r,
                                                       int protonationState,
                                                       const bead_spec_ptr_t& spec)
    {  
        auto index = indexGenerator_(*this);
        prot_bead_ptr_t bead = 
            ProtonatableBead::create(index, name, protonationState, spec);
        bead->position(r);
        protonatableBeads_.push_back(bead);
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
        return protonatableBeads_.size();
    }
                
    std::size_t CoarseGrained::protonationState() const
    {
        std::size_t numberOfProtons{0};
        for (auto pb: protonatableBeads_) {
            numberOfProtons += pb->protonationState();
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
            std::size_t index = counter + 1;
            stream.read(charBuffer, bufferSize);
            std::string name = std::string(charBuffer, bufferSize);
            boost::trim(name);
            stream.read(charBuffer, bufferSize);
            std::string specName = std::string(charBuffer, bufferSize);
            boost::trim(specName);
            bead_spec_ptr_t spec = catalog->lookup(specName);
            std::size_t id;
            stream >> id;
            if ( id != index) {
                std::string msg = 
                    "Indexing in input particle model may be wrong. Expected " + 
                    boost::lexical_cast<std::string, std::size_t>(index) + ", " +
                    "got instead " + boost::lexical_cast<std::string, std::size_t>(id) + ".";
                throw std::domain_error(msg);
            }
            real_t x, y, z, px, py, pz;
            stream >> x >> y >> z >> px >> py >> pz;
            position_t r{x, y, z};
            momentum_t p{px, py, pz};
            bool protonatable;
            stream >> protonatable;
            if ( protonatable ) {
                std::size_t protonationState;            
                stream >> protonationState;
                bead_ptr_t bead = 
                    cg->addProtonatableBead(name, r, protonationState, spec);
                bead->momentum(p);
            } else {
                bead_ptr_t bead = cg->addBead(name, r, spec);
                bead->momentum(p);
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

