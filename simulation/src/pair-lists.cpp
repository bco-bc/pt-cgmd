/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 21, 2019, 3:43 PM
 */

#include "simploce/simulation/pair-lists.hpp"
#include <utility>
#include <vector>

namespace simploce {

    PairLists::PairLists() :
        numberOfParticles_{0}, ppPairs_{}, updated_{false} {
    }

    PairLists::PairLists(std::size_t numberOfParticles,
                            std::vector<pp_pair_t> ppPairs) :
            numberOfParticles_{numberOfParticles}, ppPairs_{std::move(ppPairs)}, updated_{true} {
    }

    const PairLists::pp_pair_cont_t&
    PairLists::particlePairList() const {
        return ppPairs_;
    }

    bool
    PairLists::isModified() const {
        return updated_;
    }

    void
    PairLists::modified_(bool modified) {
        updated_ = modified;
    }

    bool
    PairLists::isEmpty() const {
        return ppPairs_.empty();
    }

    std::size_t
    PairLists::numberOfParticles() const {
        return numberOfParticles_;
    }

    void
    PairLists::numberOfParticles(std::size_t numberOfParticles) {
        numberOfParticles_= numberOfParticles;
    }

}

