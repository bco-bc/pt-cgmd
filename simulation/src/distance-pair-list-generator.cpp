/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#include "simploce/simulation/distance-pair-list-generator.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <utility>

namespace simploce {

    /**
     * Returns particle pair lists for -any- collection of particles.
     * @param box Simulation box.
     * @param cutoff Cutoff distance.
     * @param bc Boundary condition.
     * @return Particle pairs.
     */
    static typename PairLists::pp_pair_cont_t
    forParticles_(const box_ptr_t& box,
                  const dist_t cutoff,
                  const bc_ptr_t& bc,
                  const std::vector<p_ptr_t>& particles)
    {
        using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;

        static real_t rc2 = cutoff() * cutoff();

        if ( particles.empty() ) {
            return std::move(pp_pair_cont_t{});  // Empty pair list.
        }

        // Particle/particle pair list.
        pp_pair_cont_t particlePairs{};

        // For all pairs in the set.
        for (auto iter_i = particles.begin(); iter_i != (particles.end() - 1); ++iter_i) {
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
                    particlePairs.emplace_back(pair);
                }
            }
        }

        // Done.
        return std::move(particlePairs);
    }

    /**
     * Returns particle pair list for particles in -different- groups.
     * @tparam P Particle typeName.
     * @param box Simulation box.
     * @param cutoff Cutoff distance.
     * @param bc Boundary condition.
     * @param groups
     * @return Particle pairs.
     */
    static typename PairLists::pp_pair_cont_t
    forParticlesBetweenGroups_(const box_ptr_t& box,
                               const dist_t& cutoff,
                               const bc_ptr_t& bc,
                               const std::vector<pg_ptr_t>& groups)
    {
        using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
        using pp_pair_t = typename PairLists::pp_pair_t;

        static real_t rc2 =  cutoff() * cutoff();

        if ( groups.empty() ) {
            return std::move(pp_pair_cont_t{});  // Empty list.
        }

        // Pair list.
        pp_pair_cont_t particlePairs{};

        // For all different particle group pairs. Include pairs if groups are within the cutoff distance.
        for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
            const auto gi = *iter_i;
            auto r_gi = gi->position();
            const auto particles_i = gi->particles();
            for (auto iter_j = iter_i + 1; iter_j != groups.end(); ++iter_j) {
                auto gj = *iter_j;
                auto r_gj = gj->position();
                auto R = bc->apply(r_gi, r_gj);  // Group distance (vector).
                auto R2 = norm_square<real_t>(R);
                if ( R2 <= rc2 ) {
                    // Include all particle pairs between the two groups.
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

    /**
     * Returns particle pair list between any collection of particles and particles in groups.
     * @param box Simulation box.
     * @param cutoff Cutoff distance.
     * @param bc Boundary condition.
     * @param particles Particles.
     * @param groups Particle groups.
     * @return Particle pairs.
     */
    static typename PairLists::pp_pair_cont_t
    forParticlesAndGroups_(const box_ptr_t& box,
                           const dist_t cutoff,
                           const bc_ptr_t& bc,
                           const std::vector<p_ptr_t>& particles,
                           const std::vector<pg_ptr_t>& groups)
    {
        using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
        using pp_pair_t = typename PairLists::pp_pair_t;

        static real_t rc2 =  cutoff() * cutoff();

        if ( particles.empty() || groups.empty() ) {
            return std::move(pp_pair_cont_t{});  // Empty list.
        }

        // Pair list.
        pp_pair_cont_t particlePairs{};

        for (const auto& pi : particles) {
            position_t ri = pi->position();
            for (const auto& g : groups) {
                if ( !g->contains(pi) ) {
                    position_t r_g = g->position();
                    auto R = bc->apply(ri, r_g);  // Distance (vector) between particle and group.
                    auto R2 = norm_square<real_t>(R);
                    if (R2 <= rc2) {
                        for (const auto& pj : g->particles()) {
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
     * Returns particle pair list for non-bonded particle pairs for particles in the -same- group.
     * @return Particle pairs.
     */
    static typename PairLists::pp_pair_cont_t
    forParticlesInGroups(const box_ptr_t& box,
                         const bc_ptr_t& bc,
                         const std::vector<pg_ptr_t>& groups) {
        using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
        using pp_pair_t = typename PairLists::pp_pair_t;

        if ( groups.empty() ) {
            return std::move(pp_pair_cont_t{});  // Empty list.
        }

        // Pair list.
        pp_pair_cont_t particlePairs{};

        for (const auto& g: groups) {
            auto pairs = g->nonBondedParticlePairs();
            for (const auto& pair: pairs ) {
                pp_pair_t pp = std::make_pair(pair.first, pair.second);
                particlePairs.emplace_back(pp);
            }
        }

        return std::move(particlePairs);
    }

    /**
     * Make the particle pair lists.
     * @tparam P Particle typeName.
     * @param box Simulation box.
     * @param bc Boundary condition.
     * @param all All particles.
     * @param free Free particles.
     * @param groups Particle groups.
     * @return Particle pair lists.
     */
    static PairLists
    makePairLists_(const box_ptr_t& box,
                   const dist_t cutoff,
                   const bc_ptr_t& bc,
                   const std::vector<p_ptr_t> &all,
                   const std::vector<p_ptr_t> &free,
                   const std::vector<pg_ptr_t> &groups)
    {
        static util::Logger logger("simploce::makePairLists_()");
        static bool firstTime = true;
        if ( firstTime ) {
            logger.info("Creating distance-based particle pair lists.");
            logger.debug("Cutoff distance: " + util::toString(cutoff()));
        }

        // Prepare new particle pair list.
        auto particlePairs = forParticles_(box, cutoff, bc, free);
        auto ppSize = particlePairs.size();
        auto fgParticlePairs = forParticlesAndGroups_(box, cutoff, bc, free, groups);
        auto fgSize = fgParticlePairs.size();
        particlePairs.insert(particlePairs.end(), fgParticlePairs.begin(), fgParticlePairs.end());
        auto ggParticlePairs = forParticlesBetweenGroups_(box, cutoff, bc, groups);
        auto ggSize = ggParticlePairs.size();
        particlePairs.insert(particlePairs.end(), ggParticlePairs.begin(), ggParticlePairs.end());
        auto inGroupsParticlePairs = forParticlesInGroups(box, bc, groups);
        auto ingSize = inGroupsParticlePairs.size();
        particlePairs.insert(particlePairs.end(), inGroupsParticlePairs.begin(), inGroupsParticlePairs.end());

        if ( firstTime ) {
            logger.debug("Non-bonded particle pair list:");
            logger.debug("Number of free-particle/free-particle pairs: " + util::toString(ppSize));
            logger.debug("Number of free-particle/particle-in-group pairs: " + util::toString(fgSize));
            logger.debug("Number of particle-in-group/particle-in-other-group pairs: " + util::toString(ggSize));
            logger.debug("Number of particle-in-group/particle-in-same-group: " + util::toString(ingSize));
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

    DistancePairListGenerator::DistancePairListGenerator(dist_t cutoff, bc_ptr_t bc) :
        pair_lists_generator{}, cutoff_{cutoff}, bc_{std::move(bc)} {
    }

    PairLists
    DistancePairListGenerator::generate(const p_system_ptr_t& particleSystem) const {
        auto box = particleSystem->box();
        return particleSystem->doWithAllFreeGroups<PairLists>([this, box] (
                const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) {
            return std::move(makePairLists_(box, this->cutoff_, this->bc_, all, free, groups));
        });
    }

}
