/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11 October 2019, 16:29
 */

#include "simploce/simulation/cell-pair-list-generator.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/s-util.hpp"
#include "simploce/simulation/grid.hpp"
#include "simploce/simulation/cell.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include <utility>

namespace simploce {
    namespace pairlist {

        /**
         * Returns all particles that are involved in non-bonded interactions, that is
         * can form particle pairs for non-bonded interactions.
         * IMPORTANT: DOES NOT DISTINGUISH BETWEEN PARTICLE PAIRS INVOLVED IN A BONDED INTERACTION IN THE SAME
         * PARTICLE GROUP FROM THOSE PAIRS LOCATED IN DIFFERENT GROUPS. MUST BE CORRECTED.
         * @param all All particles.
         * @param free Free particles.
         * @param groups Particle groups.
         * @return Particles involved in non-bonded interactions.
         */
        static std::set<p_ptr_t> nonBondedParticles(const std::vector<p_ptr_t> &all,
                                                    const std::vector<p_ptr_t> &free,
                                                    const std::vector<pg_ptr_t> &groups) {
            util::Logger logger("simploce::pairlist::nonBondedParticles()");
            logger.trace("Entering");

            // Non-bonded particles.
            std::set<p_ptr_t> particles{};

            logger.debug(std::to_string(groups.size()) + ": Number of particle groups.");
            for (const auto &g: groups) {
                auto pairs = g->nonBondedParticlePairs();
                logger.debug(std::to_string(pairs.size()) +
                                     ": Number of particle pairs in group #" + std::to_string(g->id()));
                for (const auto &pair: pairs) {
                   particles.insert(pair.first);
                   particles.insert(pair.second);
                }
            }
            auto sizeInGroups = particles.size();
            logger.debug(std::to_string(sizeInGroups) + ": Number of non-bonded particles added from within groups.");

            // Find all particle pairs for particles in different groups.
            if (!groups.empty()) {
                for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
                    const auto &gi = *iter_i;
                    const auto particles_i = gi->particles();
                    for (const auto &pi: particles_i) {
                        particles.insert(pi);
                        for (auto iter_j = iter_i + 1; iter_j != groups.end(); ++iter_j) {
                            auto &gj = *iter_j;
                            const auto particles_j = gj->particles();
                            for (const auto &pj: particles_j) {
                                particles.insert(pj);
                            }
                        }
                    }
                }
            }
            auto sizeBetweenGroups = particles.size() - sizeInGroups;
            logger.debug(std::to_string(sizeBetweenGroups) + ": Number of non-bonded particles added from between groups.");

            // Find all particle pairs for free particles and particles in groups.
            for (const auto& pi : free) {
                for (const auto &g: groups) {
                    if (!g->contains(pi)) {
                        auto& particles_j = g->particles();
                        for (auto& pj: particles_j) {
                            particles.insert(pj);
                        }
                    }
                }
            }
            auto sizeFreeGroups = particles.size() - sizeBetweenGroups - sizeInGroups;
            logger.debug(std::to_string(sizeFreeGroups) +
                         ": Number of non-bonded particles between free particles and particles in groups.");

            // Free particles can always form pairs for non-bonded interactions.
            for (const auto& p : free) {
                particles.insert(p);
            }
            logger.debug(std::to_string(free.size()) +
                         ": Number of non-bonded particles from free particles.");
            logger.info(std::to_string(particles.size()) +
                        ": Total number of particles involved in non-bonded interactions.");

            logger.trace("Leaving.");
            return std::move(particles);
        }

        static PairLists
        makePairListNonConcurrent(const box_ptr_t &box,
                                  const dist_t & cutoff,
                                  const bc_ptr_t &bc,
                                  const std::vector<p_ptr_t>& particles,
                                  const Grid& grid) {
            using pp_pair_cont_t = typename PairLists::pp_pair_cont_t;
            using pp_pair_t = typename PairLists::pp_pair_t;

            static util::Logger logger("simploce::pairlist::makePairListNonConcurrent()");
            logger.trace("Entering.");

            static real_t rc2 = cutoff() * cutoff();
            logger.debug(std::to_string(cutoff()) + ": Cutoff distance for cell-based particle pairs.");

            // Find unique particle pairs.
            std::map<std::string, pp_pair_t> pairs{};
            auto insertParticlePair = [&pairs] (const p_ptr_t& pi, const p_ptr_t& pj) {
                // Each particle pair should appear just once in the pair list.
                auto index_i = pi->index();
                auto index_j = pj->index();
                auto key = index_j > index_i ?
                       std::to_string(index_i) + "-" + std::to_string(index_j) :
                       std::to_string(index_j) + "-" + std::to_string(index_i);
                pp_pair_t ppPair = std::make_pair(pi, pj);
                                    auto pair = std::make_pair(key, ppPair);
                                    pairs.insert(pair);
            };

            auto neighboring = grid.neighbors();
            for (auto& c : neighboring) {
                // Get central cell.
                auto location = c.first;
                auto central = grid.findCell(location);
                auto cParticles = central->particles();

                if ( !cParticles.empty() ) {
                    // Particle pairs in central cell.
                    for (auto i = 0; i < (cParticles.size() - 1); ++i) {
                        auto &pi = cParticles[i];
                        auto ri = pi->position();
                        auto index_i = pi->index();
                        for (auto j = i + 1; j < cParticles.size(); ++j) {
                            auto &pj = cParticles[j];
                            auto rj = pj->position();
                            auto rij = bc->apply(ri, rj);
                            auto rij2 = norm_square<real_t>(rij);
                            if (rij2 <= rc2) {
                                insertParticlePair(pi, pj);
                            }
                        }
                    }

                    // Particle pairs between particles in central cell and in neighboring cells.
                    auto neighbors = c.second;
                    for (const auto &pi: cParticles) {
                        auto index_i = pi->index();
                        auto ri = pi->position();
                        for (const auto &ce: neighbors) {
                            for (const auto &pj: ce->particles()) {
                                auto rj = pj->position();
                                auto rij = bc->apply(ri, rj);
                                auto rij2 = norm_square<real_t>(rij);
                                if (rij2 <= rc2) {
                                    insertParticlePair(pi, pj);
                                }
                            }
                        }
                    }
                }

                // Done.
            }

            // Pair list.
            pp_pair_cont_t particlePairs{};
            for (auto& p: pairs) {
                auto pair = p.second;
                particlePairs.emplace_back(pair);
            }

            logger.debug(std::to_string(particlePairs.size()) + ": Number of particle pairs.");
            auto total = particles.size() * (particles.size() - 1) / 2;
            logger.debug("Total number of POSSIBLE particle pairs: " + util::toString(total));
            logger.debug("Fraction (%) of total number of possible particle pairs: " +
                         util::toString(real_t(particlePairs.size()) * 100.0 / real_t(total)));

            logger.trace("Leaving.");
            return std::move(PairLists{particles.size(), particlePairs});
        }

        static PairLists
        makePairListConcurrent(const box_ptr_t &box,
                                  const dist_t & cutoff,
                                  const bc_ptr_t &bc,
                                  const std::vector<p_ptr_t>& particles,
                                  const Grid& grid) {
            static util::Logger logger("simploce::pairlist::makePairListConcurrent()");
            logger.trace("Entering.");

            auto pairLists =
                    std::move(pairlist::makePairListNonConcurrent(box, cutoff, bc, particles, grid));

            logger.trace("Leaving.");
            return std::move(pairLists);
        }

    }



    CellPairListGenerator::CellPairListGenerator(param_ptr_t param, bc_ptr_t bc) :
        param_{std::move(param)}, bc_{std::move(bc)} {
    }
    
    PairLists
    CellPairListGenerator::generate(const p_system_ptr_t& particleSystem) const
    {
        static util::Logger logger("simploce::CellPairListGenerator::generate()");
        logger.trace(("Entering."));

        static int counter = 1;
        if (counter == 1) {
            logger.info("Cell-based particle pair list generation.");
        }

        // Setup.
        static auto particles = particleSystem->doWithAllFreeGroups<std::vector<p_ptr_t>>([] (
                const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) {
            auto nbParticles = std::move(pairlist::nonBondedParticles(all, free, groups));
            std::vector<p_ptr_t> particles{};
            for (const auto& p : nbParticles) {
                particles.emplace_back(p);
            }
            return std::move(particles);
        });
        static dist_t cutoff = util::computePairListCutoff(this->param_, particleSystem);
        static Grid grid(particleSystem->box(), bc_, cutoff);

        // Generate pair list.
        grid.assignParticlesToCells(bc_, particles);
        auto pairLists =
                pairlist::makePairListConcurrent(particleSystem->box(),
                                                 cutoff,
                                                 bc_,
                                                 particles,
                                                 grid);
        if (counter == 1) {
            logger.info(std::to_string(pairLists.particlePairList().size()) +
                                ": Number of particle pairs in step #" + std::to_string(counter));
        }
        counter += 1;

        logger.trace("Leaving.");
        return std::move(pairLists);
    }

}