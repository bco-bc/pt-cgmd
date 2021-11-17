/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 21, 2019, 3:43 PM
 */

#ifndef PAIR_LISTS_HPP
#define PAIR_LISTS_HPP

#include "simploce/particle/particle-group.hpp"

namespace simploce {
    
    /**
     * Holds particle/particle, particle/particle group, and
     * particle group/particle group pair lists for non-bonded interactions.
     */
    class PairLists {
    public:
        
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
                  pp_pair_cont_t pairList);
        
        /**
         * Returns particle/particle pairs.
         * @return Particle pairs.
         */
        const pp_pair_cont_t & particlePairList() const;
               
        /**
         * Were these pair lists altered.
         * @return Result.
         */
        bool isModified() const;

        /**
         * Are there any pair lists.
         * @return Result.
         */
        bool isEmpty() const;

        /**
         * Total number of particles in the particle system.
         * @return Number.
         */
        std::size_t numberOfParticles() const;
        
    private:

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
    
}

#endif /* PAIR_LISTS_HPP */

