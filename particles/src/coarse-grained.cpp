/*
 * File:   coarse-grained.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:38 PM
 */

#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/conf/p-conf.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/particle/bond.hpp"
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
        ParticleSystem{} {
    }
    
    p_ptr_t
    CoarseGrained::addBead(const std::string& name,
                           const spec_ptr_t& spec) {
        return this->addParticle(name, spec);
    }

    std::size_t
    CoarseGrained::numberOfBeads() const {
        return this->numberOfParticles();
    }

    p_ptr_t CoarseGrained::createParticle_(const id_t& id,
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

