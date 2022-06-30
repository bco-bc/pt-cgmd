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
#include "simploce/units/units-mu.hpp"
#include <utility>

namespace simploce {
    namespace pairlist {

        static dist_t
        computeCutoff_(const param_ptr_t& param,
                       const p_system_ptr_t & particleSystem) {
            util::Logger logger("simploce::pairlist::computeCutoff()");
            logger.trace("Entering.");

            temperature_t temperature = param->get<real_t>("simulation.temperature");
            auto isMesoscale = param->get<bool>("simulation.mesoscale", false);
            stime_t dt = param->get<real_t>("simulation.timestep", 0.001);
            auto cutoff = param->get<real_t>("simulation.forces.cutoff");
            auto nPairLists = param->get<int>("simulation.npairlists");

            // Compute typical displacement.
            auto numberOfParticles = particleSystem->numberOfParticles();
            auto KB = isMesoscale ? 1 : units::mu<real_t>::KB;
            mass_t mass = particleSystem->mass();
            mass /= real_t(numberOfParticles);
            energy_t eKin = KB * temperature() / mass();  // Equipartition, 1D.
            auto v = std::sqrt(2.0 * eKin() / mass());
            auto displacement = v * dt;

            // Verlet cutoff distance for pair lists.
            // See https://en.wikipedia.org/wiki/Verlet_list
            dist_t pairListCutoff = cutoff + 2.0 * nPairLists * displacement();

            logger.debug(util::toString(temperature) + ": Temperature.");
            logger.debug(util::toString(dt) + ": Time step.");
            logger.debug(util::toString(eKin) +
                          ": Estimated kinetic energy per particle at the given temperature.");
            logger.debug(util::toString(mass) + ": Average mass.");
            logger.debug(std::to_string(isMesoscale) + ": Mesoscale simulation?");
            logger.debug(std::to_string(cutoff) + ": Cutoff distance for forces.");
            logger.debug(std::to_string(v) + ": Average velocity.");
            logger.debug(std::to_string(displacement()) + ": Typical displacement in single step.");
            logger.debug(util::toString(nPairLists) + ": Number of steps between update pair lists.");
            logger.info(util::toString(pairListCutoff) + ": Cutoff distance for pair lists generation.");

            logger.trace("Leaving.");
            return pairListCutoff;
        }

        /**
         * Returns particle pair lists for -any- collection of particles.
         * @param box Simulation box.
         * @param cutoff Cutoff distance.
         * @param bc Boundary condition.
         * @return Particle pairs.
         */
        static typename PairLists::pp_pair_cont_t
        forParticles_(const box_ptr_t &box,
                      const dist_t cutoff,
                      const bc_ptr_t &bc,
                      const std::vector<p_ptr_t> &particles) {
            using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;

            static real_t rc2 = cutoff() * cutoff();

            if (particles.empty()) {
                return std::move(pp_pair_cont_t{});  // Empty pair list.
            }

            // Particle/particle pair list.
            pp_pair_cont_t particlePairs{};

            // For all pairs in the set.
            for (auto iter_i = particles.begin(); iter_i != (particles.end() - 1); ++iter_i) {
                const auto &pi = *iter_i;
                position_t ri = pi->position();
                for (auto iter_j = iter_i + 1; iter_j != particles.end(); ++iter_j) {
                    const auto &pj = *iter_j;
                    position_t rj = pj->position();
                    dist_vect_t rij = bc->apply(ri, rj);
                    auto R2 = norm_square<real_t>(rij);
                    if (R2 <= rc2) {
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
        forParticlesBetweenGroups_(const box_ptr_t &box,
                                   const dist_t &cutoff,
                                   const bc_ptr_t &bc,
                                   const std::vector<pg_ptr_t> &groups) {
            using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
            using pp_pair_t = typename PairLists::pp_pair_t;

            static real_t rc2 = cutoff() * cutoff();

            if (groups.empty()) {
                return std::move(pp_pair_cont_t{});  // Empty list.
            }

            // Pair list.
            pp_pair_cont_t particlePairs{};

            for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
                const auto& gi = *iter_i;
                const auto particles_i = gi->particles();
                for (const auto& pi : particles_i) {
                    auto ri = pi->position();
                    for (auto iter_j = iter_i + 1; iter_j != groups.end(); ++iter_j) {
                        auto &gj = *iter_j;
                        const auto particles_j = gj->particles();
                        for (const auto& pj : particles_j) {
                            auto rj = pj->position();
                            auto rij = bc->apply(ri, rj);
                            auto R2 = norm_square<real_t>(rij);
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
         * Returns particle pair list between -any- collection of particles and particles in groups.
         * @param box Simulation box.
         * @param cutoff Cutoff distance.
         * @param bc Boundary condition.
         * @param particles Particles.
         * @param groups Particle groups.
         * @return Particle pairs.
         */
        static typename PairLists::pp_pair_cont_t
        forParticlesAndGroups_(const box_ptr_t &box,
                               const dist_t cutoff,
                               const bc_ptr_t &bc,
                               const std::vector<p_ptr_t> &particles,
                               const std::vector<pg_ptr_t> &groups) {
            using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
            using pp_pair_t = typename PairLists::pp_pair_t;

            static real_t rc2 = cutoff() * cutoff();

            if (particles.empty() || groups.empty()) {
                return std::move(pp_pair_cont_t{});  // Empty list.
            }

            // Pair list.
            pp_pair_cont_t particlePairs{};

            for (const auto& pi : particles) {
                position_t ri = pi->position();
                for (const auto &g: groups) {
                    if (!g->contains(pi)) {
                        auto& particles_j = g->particles();
                        for (auto& pj: particles_j) {
                            auto rj = pj->position();
                            auto rij = bc->apply(ri, rj);
                            auto R2 = norm_square<real_t>(rij);
                            if (R2 <= rc2) {
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
        forParticlesInGroups(const box_ptr_t &box,
                             const bc_ptr_t &bc,
                             const std::vector<pg_ptr_t> &groups) {
            using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
            using pp_pair_t = typename PairLists::pp_pair_t;

            if (groups.empty()) {
                return std::move(pp_pair_cont_t{});  // Empty list.
            }

            // Pair list.
            pp_pair_cont_t particlePairs{};

            // Include all non-bonded pairs.
            for (const auto &g: groups) {
                auto pairs = g->nonBondedParticlePairs();
                for (const auto &pair: pairs) {
                    pp_pair_t pp = std::make_pair(pair.first, pair.second);
                    particlePairs.emplace_back(pp);
                }
            }

            return std::move(particlePairs);
        }

        /**
         * Generates pair lists concurrently. There must be at least 4 cores/threads available.
         * @return Particle pairs.
         */
        static PairLists
        makePairListsConcurrent_(const box_ptr_t &box,
                                 const dist_t& cutoff,
                                 const bc_ptr_t &bc,
                                 const std::vector<p_ptr_t> &all,
                                 const std::vector<p_ptr_t> &free,
                                 const std::vector<pg_ptr_t> &groups) {
            using result_t = PairLists::pp_pair_cont_t;

            static util::Logger logger("simploce::pairlist::makePairListsConcurrent_()");
            logger.trace("Entering");

            std::vector<result_t> results{};
            std::vector<std::future<result_t>> futures{};

            futures.emplace_back(std::async(
                    std::launch::async,
                    forParticles_,
                    std::ref(box),
                    std::ref(cutoff),
                    std::ref(bc),
                    std::ref(free)
            ));
            futures.emplace_back(std::async(
                    std::launch::async,
                    forParticlesAndGroups_,
                    std::ref(box),
                    std::ref(cutoff),
                    std::ref(bc),
                    std::ref(free),
                    std::ref(groups)
            ));
            futures.emplace_back(std::async(
                    std::launch::async,
                    forParticlesBetweenGroups_,
                    std::ref(box),
                    std::ref(cutoff),
                    std::ref(bc),
                    std::ref(groups)
            ));
            // Wait for tasks to complete.
            results = util::waitForAll<result_t>(futures);

            // One remaining pair list to be completed in the current thread.
            const auto result = forParticlesInGroups(box, bc, groups);
            results.emplace_back(result);

            // Collect all results.
            logger.debug("Particle pair lists (concurrent):");
            result_t particlePairs{};
            for (auto &r: results) {
                logger.debug(std::to_string(r.size()) + ": Number of particle pairs.");
                particlePairs.insert(particlePairs.end(), r.begin(), r.end());
            }

            logger.debug(std::to_string(particlePairs.size()) + ": Number of particle pairs.");
            auto total = all.size() * (all.size() - 1) / 2;
            logger.debug("Total number of POSSIBLE particle pairs: " + util::toString(total));
            logger.debug("Fraction (%) of total number of possible particle pairs: " +
                         util::toString(real_t(particlePairs.size()) * 100.0 / real_t(total)));

            logger.trace("Leaving.");
            return std::move(PairLists(all.size(), particlePairs));
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
        makePairListsNonConcurrent_(const box_ptr_t &box,
                                    const dist_t & cutoff,
                                    const bc_ptr_t &bc,
                                    const std::vector<p_ptr_t> &all,
                                    const std::vector<p_ptr_t> &free,
                                    const std::vector<pg_ptr_t> &groups) {
            static util::Logger logger("simploce::pairlist::makePairListsNonConcurrent_()");
            logger.trace("Entering.");

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

            // Some debugging information.
            logger.debug("Particle pair lists (nonconcurrent):");
            logger.debug("Cutoff distance: " + util::toString(cutoff()));
            logger.debug("Number of free-particle/free-particle pairs: " + util::toString(ppSize));
            logger.debug("Number of free-particle/particle-in-group pairs: " + util::toString(fgSize));
            logger.debug("Number of particle-in-group/particle-in-other-group pairs: " + util::toString(ggSize));
            logger.debug("Number of particle-in-group/particle-in-same-group: " + util::toString(ingSize));
            logger.debug("Total number of particle pairs: " + util::toString(particlePairs.size()));
            auto total = all.size() * (all.size() - 1) / 2;
            logger.debug("Total number of POSSIBLE particle pairs: " + util::toString(total));
            logger.debug("Fraction (%) of total number of possible particle pairs: " +
                         util::toString(real_t(particlePairs.size()) * 100.0 / real_t(total)));

            logger.trace("Leaving.");
            return std::move(PairLists(all.size(), particlePairs));
        }

    }

    DistancePairListGenerator::DistancePairListGenerator(param_ptr_t param, bc_ptr_t bc) :
        pair_lists_generator{}, param_{std::move(param)}, bc_{std::move(bc)} {
    }

    PairLists
    DistancePairListGenerator::generate(const p_system_ptr_t& particleSystem) const {
        auto box = particleSystem->box();
        static dist_t cutoff = pairlist::computeCutoff_(this->param_, particleSystem);
        return particleSystem->doWithAllFreeGroups<PairLists>([this, box] (
                const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) {

            if (std::thread::hardware_concurrency() >= 4) {
                return std::move(pairlist::makePairListsConcurrent_(box, cutoff, this->bc_, all, free, groups));
            } else {
                return std::move(pairlist::makePairListsNonConcurrent_(box, cutoff, this->bc_, all, free, groups));
            }
        });
    }

}
