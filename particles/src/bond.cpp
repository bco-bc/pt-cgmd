/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on November 16, 2021.
 */

#include "simploce/particle/bond.hpp"
#include "simploce/particle/particle.hpp"
#include <stdexcept>
#include <memory>
#include <utility>

namespace simploce {

    Bond
    Bond::makeBond(const p_ptr_t& p1,
                   const p_ptr_t& p2) {
        return Bond{p1, p2};
    }

    p_ptr_t
    Bond::getParticleOne() const {
        return p1_;
    }

    p_ptr_t
    Bond::getParticleTwo() const {
        return p2_;
    }

    bool
    Bond::contains(const p_ptr_t& particle) const {
      return particle == p1_ || particle == p2_;
    }

    Bond::Bond(p_ptr_t  p1,
               p_ptr_t  p2) :
      p1_{std::move(p1)}, p2_{std::move(p2)} {
        if ( !p1_ || !p2_ ) {
            throw std::domain_error("Bond: Two particles must be provided.");
        }
        if ( p1_ == p2_ || p1_->id() == p2_->id() ) {
            throw std::domain_error("Bond: Particles must be different.");
        }
    }

}

