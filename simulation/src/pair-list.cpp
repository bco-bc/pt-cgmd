/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 21, 2019, 3:43 PM
 */

#include "simploce/simulation/pair-list.hpp"
#include <utility>
#include <vector>

namespace simploce {

    using p_pair_t = PairList::p_pair_t;

    PairList::PairList() :
        numberOfParticles_{0}, nbParticlePairs_{}, bParticlePairs_{}, modified_{false} {
    }

    const std::vector<p_pair_t>&
    PairList::nonBondedParticlePairs() const {
        return nbParticlePairs_;
    }

    const std::vector<p_pair_t>&
    PairList::bondedParticlePairs() const {
        return bParticlePairs_;
    }

    bool
    PairList::isModified() const {
        return modified_;
    }

    bool
    PairList::isEmpty() const {
        return nbParticlePairs_.empty();
    }

    std::size_t
    PairList::numberOfParticles() const {
        return numberOfParticles_;
    }

    void
    PairList::modified(bool modified) {
        modified_ = modified;
    }

    void
    PairList::numberOfParticles(std::size_t numberOfParticles) {
        numberOfParticles_= numberOfParticles;
    }

    void
    PairList::nonBoundedParticlePairs(const std::vector<p_pair_t> &particlePairs) {
        nbParticlePairs_ = particlePairs;
    }

    void
    PairList::bondedParticlePairs(const std::vector<p_pair_t> &particlePairs) {
        bParticlePairs_ = particlePairs;
    }

}

