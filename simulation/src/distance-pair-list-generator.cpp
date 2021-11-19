/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#include "simploce/simulation/distance-pair-list-generator.hpp"

#include <utility>
#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"

namespace simploce {

    /**
     * For -any- collection of particles, compute the pair lists.
     * @param box Simulation box.
     * @param bc Boundary condition.
     * @return Particle pairs.
     */
    static typename PairLists::pp_pair_cont_t
    forParticles_(const box_ptr_t& box,
                  const bc_ptr_t& bc,
                  const std::vector<p_ptr_t>& particles)
    {
        using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;

        static real_t rc2 = properties::squareCutoffDistance(box);

        if ( particles.empty() ) {
            return pp_pair_cont_t{};  // Empty pair list.
        }

        // Particle/particle pair list.
        pp_pair_cont_t particlePairs{};

        // For all pairs in the set.
        for (auto iter_i = particles.begin(); iter_i != particles.end() - 1; ++iter_i) {
            const auto& pi = *iter_i;
            position_t ri = pi->position();
            for (auto iter_j = iter_i + 1; iter_j != particles.end(); ++iter_j) {
                const auto& pj = *iter_j;
                position_t rj = pj->position();
                dist_vect_t R = bc->apply(ri, rj);
                auto R2 = norm_square<real_t>(R);
                if ( R2 <= rc2 ) {
                    // Include this pair.
                    auto pair = std::make_pair(pi, pj);
                    particlePairs.push_back(pair);
                }
            }
        }

        // Done.
        return std::move(particlePairs);
    }

    /**
     * Particle group pair list.
     * @tparam P Particle type.
     * @param box Simulation box.
     * @param bc Boundary condition.
     * @param groups
     * @return Particle group pairs.
     */
    static typename PairLists::pp_pair_cont_t
    forGroups_(const box_ptr_t& box,
               const bc_ptr_t& bc,
               const std::vector<pg_ptr_t>& groups)
    {
        using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
        using pp_pair_t = typename PairLists::pp_pair_t;

        static real_t rc2 =  properties::squareCutoffDistance(box);

        if ( groups.empty() ) {
            return pp_pair_cont_t{};  // Empty list.
        }

        // Pair list.
        pp_pair_cont_t particlePairs{};

        // For all particle group pairs.
        for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
            const auto gi = *iter_i;
            auto r_gi = gi->position();
            const auto particles_i = gi->particles();
            for (auto iter_j = iter_i + 1; iter_j != groups.end(); ++iter_j) {
                auto gj = *iter_j;
                auto r_gj = gj->position();
                auto R = bc->apply(r_gi, r_gj);
                auto R2 = norm_square<real_t>(R);
                if ( R2 <= rc2 ) {
                    // Include all particle pairs.
                    const auto particles_j = gj->particles();
                    for (const auto& pi : particles_i) {
                        for (const auto& pj : particles_j) {
                            pp_pair_t pair = std::make_pair(pi, pj);
                            particlePairs.push_back(pair);
                        }
                    }
                }
            }
        }

        // Done.
        return std::move(particlePairs);
    }

    static typename PairLists::pp_pair_cont_t
    forParticlesAndGroups_(const box_ptr_t& box,
                           const bc_ptr_t& bc,
                           const std::vector<p_ptr_t>& particles,
                           const std::vector<pg_ptr_t>& groups)
    {
        using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
        using pp_pair_t = typename PairLists::pp_pair_t;

        static real_t rc2 =  properties::squareCutoffDistance(box);

        if ( particles.empty() || groups.empty() ) {
            return pp_pair_cont_t{};  // Empty list.
        }

        // Pair list.
        pp_pair_cont_t particlePairs{};

        for (const auto& pi : particles) {
            position_t ri = pi->position();
            for (const auto& g : groups) {
                if ( !g->contains(pi) ) {
                    for (const auto& pj : g->particles()) {
                        auto rj = pj->position();
                        auto R = bc->apply(ri, rj);
                        auto R2 = norm_square<real_t>(R);
                        if ( R2 <= rc2 ) {
                            pp_pair_t pair = std::make_pair(pi, pj);
                            particlePairs.push_back(pair);
                        }
                    }
                }
            }
        }

        // Done.
        return std::move(particlePairs);
    }

    /**
     * Make the particle pair lists.
     * @tparam P Particle type.
     * @param box Simulation box.
     * @param bc Boundary condition.
     * @param all All particles.
     * @param free Free particles.
     * @param groups Particle groups.
     * @return Particle pair lists.
     */
    static PairLists
    makePairLists_(const box_ptr_t& box,
                   const bc_ptr_t& bc,
                   const std::vector<p_ptr_t> &all,
                   const std::vector<p_ptr_t> &free,
                   const std::vector<pg_ptr_t> &groups)
    {
        static util::Logger logger("simploce::makePairLists_()");
        static bool firstTime = true;
        if ( firstTime ) {
            logger.info("Creating distance-based particle pair lists.");
            logger.debug("Cutoff distance: " + util::toString(properties::cutoffDistance(box)));
        }

        // Prepare new particle pair list.
        auto particlePairs = forParticles_(box, bc, free);
        auto ppSize = particlePairs.size();
        auto fgParticlePairs = forParticlesAndGroups_(box, bc, free, groups);
        auto fgSize = fgParticlePairs.size();
        particlePairs.insert(particlePairs.end(), fgParticlePairs.begin(), fgParticlePairs.end());
        auto ggParticlePairs = forGroups_(box, bc, groups);
        auto ggSize = ggParticlePairs.size();
        particlePairs.insert(particlePairs.end(), ggParticlePairs.begin(), ggParticlePairs.end());

        if ( firstTime ) {
            logger.debug("Number of free-particle/free-particle pairs: " + util::toString(ppSize));
            logger.debug("Number of free-particle/particle-in-group pairs: " + util::toString(fgSize));
            logger.debug("Number of particle-in-group/particle-in-group pairs: " + util::toString(ggSize));
            logger.debug("Total number of particle pairs: " + util::toString(particlePairs.size()));
            auto total = all.size() * (all.size() - 1) / 2;
            logger.debug("Total number of POSSIBLE particle pairs: " + util::toString(total));
            logger.debug("Fraction (%) of total number of possible particle pairs: " +
                          util::toString(real_t(particlePairs.size()) * 100.0 / real_t(total)));
            firstTime = false;
        }

        // Done.
        return std::move(PairLists(all.size(), particlePairs));
    }

    DistancePairListGenerator::DistancePairListGenerator(bc_ptr_t bc) :
        pair_lists_generator{}, bc_{std::move(bc)} {
    }

    PairLists
    DistancePairListGenerator::generate(const p_system_ptr_t& particleSystem) const {
        auto box = particleSystem->box();
        return particleSystem->doWithAllFreeGroups<PairLists>([this, box] (
                const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) {
            return std::move(makePairLists_(box, this->bc_, all, free, groups));
        });
    }

}
