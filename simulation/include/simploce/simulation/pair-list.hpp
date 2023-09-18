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
     * Holds particle pair lists for both bonded and non-bonded particles.
     */
    class PairList {
    public:
        
        /**
         * Particle pair type.
         */
        using p_pair_t = std::pair<p_ptr_t, p_ptr_t>;

        /**
         * Default constructor. Empty pair list.
         */
        PairList();

        ~PairList() = default;

        /**
         * Returns non-bonded particle pairs.
         * @return Particle pairs.
         */
        const std::vector<p_pair_t>& nonBondedParticlePairs() const;

        /**
         * Returns bonded particle pairs.
         * @return Particle pairs.
         */
        const std::vector<p_pair_t>& bondedParticlePairs() const;
               
        /**
         * Was this pair list altered?
         * @return Result.
         */
        bool isModified() const;

        /**
         * Is this list empty?
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
         * Specifies whether this pair list was modified.
         * @param modified True, if the list was updated, otherwise false.
         */
        void modified(bool modified);

        /**
         * Sets total number of particles in particle system.
         * @param numberOfParticles Number of particles.
         */
        void numberOfParticles(std::size_t numberOfParticles);

        /**
         * Sets/updates non-bonded particle pairs.
         * @param particlePairs Particle pairs.
         */
        void nonBoundedParticlePairs(const std::vector<p_pair_t>& particlePairs);

        /**
         * Sets/updates bonded particle pairs.
         * @param particlePairs Particle pairs.
         */
        void bondedParticlePairs(const std::vector<p_pair_t>& particlePairs);

        std::size_t numberOfParticles_;
        std::vector<p_pair_t> nbParticlePairs_;  // non-bonded
        std::vector<p_pair_t> bParticlePairs_;   // bonded
        bool modified_;
    };
    
}

#endif /* PAIR_LISTS_HPP */

