/*
 * File:   atomistic.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:27 PM
 */

#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/p-factory.hpp"

namespace simploce {

    at_sys_ptr_t
    Atomistic::obtainFrom(std::istream& stream, const spec_catalog_ptr_t& catalog) {
        auto atomistic = factory::atomistic();
        atomistic->parse(stream, catalog);
        return atomistic;
    }

    Atomistic::Atomistic() : ParticleSystem{} {
    }
    
    p_ptr_t
    Atomistic::addAtom(const std::string& name,
                       const spec_ptr_t& spec)
    {
        return this->addParticle(name, spec);
    }
    
    std::size_t
    Atomistic::numberOfAtoms() const
    {
        return this->numberOfParticles();
    }

    p_ptr_t
    Atomistic::createParticle_(const id_t& id,
                               int index,
                               const std::string& name,
                               const spec_ptr_t& spec) {
        return Atom::create(id, index, name, spec);
    }

    std::ostream&
    operator << (std::ostream& stream, const Atomistic& atomistic) {
        atomistic.write(stream);
        return stream;
    }
    
}

