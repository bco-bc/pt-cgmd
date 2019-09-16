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
 * File:   atomistic.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:27 PM
 */

#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/protonation-site.hpp"
#include "simploce/particle/ptypes.hpp"
#include "simploce/particle/protonation-site-catalog.hpp"
#include <memory>

namespace simploce {
    
    Atomistic::Atomistic() : ParticleModel<Atom,atom_group_t>{}, protonationSites_{}
    {        
    }
    
    atom_ptr_t Atomistic::addAtom(const std::string& name,
                                  const position_t& r, 
                                  const atom_spec_ptr_t& spec)
    {
        int index = this->size();
        atom_ptr_t atom = Atom::create(index, name, spec);
        atom->position(r);
        this->addFree(atom);
        return atom;
    }
    
    void Atomistic::protonationSites(const prot_site_catalog_ptr_t& catalog)
    {
        protonationSites_ = catalog->lookup(*this);
    }
    
    std::size_t Atomistic::numberOfProtonationSites() const
    {
        return protonationSites_.size();
    }
    
    
    std::size_t Atomistic::protonationState() const
    {
        std::size_t numberOfProtons{0};
        for (auto ps : protonationSites_) {
            numberOfProtons += ps->protonationState();
        }
        return numberOfProtons;
    }    
    
}

