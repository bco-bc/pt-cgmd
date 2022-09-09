/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#include "simploce/simulation/distance-pair-list-generator.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-util.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <utility>
#include <thread>

namespace simploce {
    namespace pairlist {

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

            static util::Logger logger("simploce::pairlist::forParticles_()");
            logger.trace("Entering.");

            static real_t rc2 = cutoff() * cutoff();

            // Empty particle set?
            if (particles.empty()) {
                // Return empty pair list.
                return std::move(pp_pair_cont_t{});
            }

            // Pair list.
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
            logger.trace("Leaving.");
            return std::move(particlePairs);
        }

        /**
         * Finds all pairs between two sets of particles.
         * @param box Simulation box
         * @param cutoff Cutoff distance.
         * @param bc Boundary condition.
         * @param particles1 Particle set #1.
         * @param particles2 Particle set #2. May be identical to particle set #1.
         * @return Particle pairs.
         */
        static typename PairLists::pp_pair_cont_t
        forParticlesParticles_(const box_ptr_t &box,
                               const dist_t &cutoff,
                               const bc_ptr_t &bc,
                               const std::vector<p_ptr_t>& particles1,
                               const std::vector<p_ptr_t>& particles2) {
            using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;

            static util::Logger logger("simploce::pairlist::forParticlesParticles_()");
            logger.trace("Entering.");

            static int counter = 1;
            counter += 1;
            std::cout << "COUNTER: " << counter << std::endl;

            logger.debug(std::to_string(particles1.size()) + ": Number of particles in particle set 1.");
            logger.debug(std::to_string(particles2.size()) + ": Number of particles in particle set 2.");
            logger.debug(std::to_string(particles1 == particles2) + ": Identical particles sets?");

            // Empty particle set(s)?
            if (particles1.empty() || particles2.empty()) {
                // Return an empty pair list.
                return std::move(pp_pair_cont_t{});
            }

            // Same set of particles?
            if ( &particles1 == &particles2) {
                return forParticles_(box, cutoff, bc, particles1);
            }

            static real_t rc2 = cutoff() * cutoff();

            pp_pair_cont_t particlePairs{};
            for (auto& pi : particles1) {
                auto ri = pi->position();
                for (auto& pj : particles2) {
                    auto rj = pj->position();
                    dist_vect_t rij = bc->apply(ri, rj);
                    auto R2 = norm_square<real_t>(rij);
                    if (R2 <= rc2) {
                        // Include this pair.
                        auto pair = std::make_pair(pi, pj);
                        particlePairs.emplace_back(pair);
                    }
                }

            }

            logger.trace("Leaving.");
            return std::move(particlePairs);
        }


        /**
         * Generates concurrently particle pair list for -any- collection of particles.
         * @param box Simulation box.
         * @param cutoff Cutoff distance.
         * @param bc Boundary condition.
         * @return Particle pairs.
         */
        static typename PairLists::pp_pair_cont_t
        forParticlesConcurrent_(const box_ptr_t &box,
                                const dist_t &cutoff,
                                const bc_ptr_t &bc,
                                const std::vector<p_ptr_t> &particles) {
            using result_t = PairLists::pp_pair_cont_t;

            static util::Logger logger("simploce::pairlist::forParticlesConcurrent_()");
            logger.trace("Entering.");

            logger.debug(std::to_string(particles.size()) + ": Number of particles.");

            // Stores results of individual tasks.
            std::vector<result_t> results{};

            // Can we use multiple threads.
            int numberOfTasks = int(std::thread::hardware_concurrency() - 4);
            numberOfTasks = numberOfTasks > 1 ? numberOfTasks : 1;
            if (numberOfTasks > 1) {
                // Yes: Generate pair list concurrently, that is, use different threads.
                int numberOfGroups = int(std::sqrt(numberOfTasks));
                int numberParticlesPerGroup = int(particles.size()/numberOfGroups);
                logger.debug(std::to_string(numberOfTasks) + ": Maximum number of tasks.");
                logger.debug(std::to_string(numberOfGroups) + ": Estimated number of groups.");
                logger.debug(std::to_string(numberParticlesPerGroup) + ": Estimated number of particles per group.");

                // Divide particles in groups.
                std::vector<std::vector<p_ptr_t>> groups{};
                std::vector<p_ptr_t> group{};
                int counter = 0;
                for (auto& p : particles) {
                    if (counter == numberParticlesPerGroup) {
                        logger.debug(std::to_string(group.size()) +
                                     ": Number of particles in current group.");
                        groups.emplace_back(group);
                        group.clear();
                        counter = 0;
                    }
                    group.emplace_back(p);
                    counter += 1;
                }
                logger.debug(std::to_string(group.size()) +
                             ": Number of particles in last group.");
                groups.emplace_back(group);
                logger.debug(std::to_string(groups.size()) + ": Actual number of groups.");

                std::vector<std::future<result_t>> futures{};
                for (auto iter_i = groups.begin(); iter_i < groups.end() - 1; ++iter_i) {
                    // Exclude last group, is handled in the current thread.
                    for (auto iter_j = iter_i; iter_j < groups.end(); ++iter_j) {
                        futures.push_back(
                                std::async(
                                        std::launch::async,
                                        forParticlesParticles_,
                                        std::ref(box),
                                        std::ref(cutoff),
                                        std::ref(bc),
                                        std::ref(*iter_i),
                                        std::ref(*iter_j)
                                )
                        );
                    }
                }
                results = util::waitForAll(futures);
                auto &particles_1 = *(groups.end() - 1);
                auto result = forParticlesParticles_(box, cutoff, bc, particles_1, particles_1);
                results.emplace_back(result);
            } else {
                // NO: Use the current thread.
                auto result = forParticles_(box, cutoff, bc, particles);
                results.emplace_back(result);
            }

            // Make the pair list.
            PairLists::pp_pair_cont_t particlePairs{};
            for (auto& r : results) {
                for (auto& pair : r) {
                    particlePairs.emplace_back(pair);
                }
            }

            logger.trace("Leaving.");
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

            // Four tasks.
            futures.emplace_back(std::async(
                    std::launch::async,
                    forParticlesConcurrent_,
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

            // One remaining task to be completed in the current thread.
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
            logger.debug("Particle pair lists (non-concurrent):");
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
        static util::Logger logger("simploce::DistancePairListGenerator::generate()");
        logger.trace("Entering.");

        static int counter = 1;
        if (counter == 1) {
            logger.info("Distance-based particle pair list generation.");
        }
        auto box = particleSystem->box();
        static dist_t cutoff = util::computePairListCutoff(this->param_, particleSystem);
        auto pairLists = particleSystem->doWithAllFreeGroups<PairLists>([this, box] (
                const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) {

            if (std::thread::hardware_concurrency() >= 4) {
                return std::move(pairlist::makePairListsConcurrent_(box, cutoff, this->bc_, all, free, groups));
            } else {
                return std::move(pairlist::makePairListsNonConcurrent_(box, cutoff, this->bc_, all, free, groups));
            }
        });
        if (counter == 1) {
            logger.info(std::to_string(pairLists.particlePairList().size()) +
                                ": Number of particle pairs in step #" + std::to_string(counter));
        }
        counter += 1;

        logger.trace("Leaving.");
        return std::move(pairLists);
    }

}
