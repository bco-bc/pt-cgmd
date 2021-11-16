/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 21, 2019, 3:43 PM
 */

#ifndef PAIR_LISTS_HPP
#define PAIR_LISTS_HPP

#include "simploce/particle/particle-group.hpp"
#include <utility>
#include <memory>
#include <vector>

namespace simploce {
    
    /**
     * Holds particle/particle, particle/particle group, and
     * particle group/particle group pair lists for non-bonded interactions.
     * @param P Particle type.
     */
    template <typename P>
    class PairLists {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Particle group pointer type.
         */
        //using pg_ptr_t = std::shared_ptr<ParticleGroup<P>>;
        
        /**
         * Particle pair type.
         */
        using pp_pair_t = std::pair<p_ptr_t, p_ptr_t>;

        /**
         * Particle pairs container type.
         */
        using pp_pair_cont_t = std::vector<pp_pair_t>;
        
        /**
         * Default constructor. Empty pair lists.
         */
        PairLists();
        
        /**
         * Constructor.
         * @param numberOfParticles Total number of particles.
         * @param Particle/particle pair list.
         */
        PairLists(std::size_t numberOfParticles,
                  const pp_pair_cont_t &pairList);
        
        /**
         * Returns particle/particle pairs.
         * @return Particle pairs.
         */
        const pp_pair_cont_t & particlePairList() const {
            return ppPairs_;
        }
               
        /**
         * Were these pair lists altered.
         * @return Result.
         */
        bool isModified() const { return updated_; }

        bool isEmpty() const { return ppPairs_.empty(); }

        std::size_t numberOfParticles() const;
        
    private:
        
        template <typename PP>
        friend class Interactor;
        
        /**
         * Signals whether the pair list was modified.
         * @param modified If true, list was updated, otherwise not.
         */
        void modified_(bool modified);

        void numberOfParticles(std::size_t numberOfParticles);

        std::size_t numberOfParticles_;
        pp_pair_cont_t ppPairs_;
        bool updated_;
    };
    
    template <typename P>
    PairLists<P>::PairLists() :
        numberOfParticles_{0}, ppPairs_{}, updated_{false} {
    }
        
    template <typename P>
    PairLists<P>::PairLists(std::size_t numberOfParticles,
                            const std::vector<pp_pair_t>& pairList) :
        numberOfParticles_{numberOfParticles}, ppPairs_{pairList}, updated_{true} {
    }
        
    template <typename P>
    void 
    PairLists<P>::modified_(bool modified) {
        updated_ = modified;
    }

    template <typename P>
    std::size_t
    PairLists<P>::numberOfParticles() const {
        return numberOfParticles_;
    }

    template <typename P>
    void
    PairLists<P>::numberOfParticles(std::size_t numberOfParticles) {
        numberOfParticles_= numberOfParticles;
    }
}

#endif /* PAIR_LISTS_HPP */

