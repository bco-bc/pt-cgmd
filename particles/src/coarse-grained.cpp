/*
 * File:   coarse-grained.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:38 PM
 */

#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/p-conf.hpp"
#include "simploce/particle/particle-group.hpp"
#include <vector>
#include <utility>
#include <memory>

namespace simploce {

    cg_sys_ptr_t
    CoarseGrained::obtainFrom(std::istream& stream,
                              const spec_catalog_ptr_t& catalog) {
        auto coarseGrained = factory::coarseGrained();
        coarseGrained->parse(stream, catalog);
        return coarseGrained;
    }

    CoarseGrained::CoarseGrained() :
            ParticleSystem<Bead,bead_group_t>{} {
    }
    
    bead_ptr_t 
    CoarseGrained::addBead(const std::string& name,
                           const spec_ptr_t& spec) {
        return this->addParticle(name, spec);
    }

    bead_group_ptr_t 
    CoarseGrained::addBeadGroup(const std::vector<bead_ptr_t>& beads, 
                                const std::vector<id_pair_t>& bonds) {
        bead_group_ptr_t group = std::make_shared<ParticleGroup<Bead>>(beads, bonds);
        this->addGroup(group);
        return group;
    }

    std::size_t 
    CoarseGrained::numberOfBeads() const {
        return this->numberOfParticles();
    }

    bead_ptr_t CoarseGrained::createParticle_(const id_t& id,
                                              int index,
                                              const std::string& name,
                                              const spec_ptr_t& spec) {
        return Bead::create(id, index, name, spec);
    }


    std::ostream&
    operator << (std::ostream& stream, const CoarseGrained& cg) {
        cg.write(stream);
        return stream;
    }
    
}

