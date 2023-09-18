/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11 October 2019, 16:29
 */

#include "simploce/simulation/cell-pair-list-generator.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-util.hpp"
#include "simploce/simulation/grid.hpp"
#include "simploce/simulation/cell.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/bond.hpp"
#include "simploce/util/logger.hpp"
#include <utility>

namespace simploce {
    namespace pairlist {

        /**
         * Returns all particles that can form particle pairs for non-bonded interactions.
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

            // Within particle groups.
            logger.debug(std::to_string(groups.size()) + ": Number of particle groups.");
            for (const auto &g: groups) {
                // Only non-bonded particles.
                auto pairs = g->nonBondedParticlePairs();
                for (const auto &pair: pairs) {
                   particles.insert(pair.first);
                   particles.insert(pair.second);
                }
            }
            auto sizeInGroups = particles.size();
            logger.debug(std::to_string(sizeInGroups) +
                         ": Number of non-bonded particles added from within groups.");

            // Find all particle pairs for particles in different groups.
            if (!groups.empty()) {
                for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
                    const auto &gi = *iter_i;
                    std::clog << "Group: " << gi->id() << std::endl;
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
            logger.debug(std::to_string(sizeBetweenGroups) +
                         ": Number of non-bonded particles added from between groups.");

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

        static std::string generatePairKey(const p_ptr_t& pi, const p_ptr_t& pj) {
            auto index_i = pi->index();
            auto index_j = pj->index();
            return index_j > index_i ?
                   std::to_string(index_i) + "-" + std::to_string(index_j) :
                   std::to_string(index_j) + "-" + std::to_string(index_i);
        }

        static std::vector<PairList::p_pair_t>
        makePairListNonConcurrent(const dist_t & cutoff,
                                  const bc_ptr_t &bc,
                                  const std::vector<p_ptr_t>& particles,
                                  const std::vector<pg_ptr_t> &groups,
                                  const Grid& grid,
                                  bool excludeFrozen) {
            static util::Logger logger("simploce::pairlist::makePairListNonConcurrent()");
            logger.trace("Entering.");

            static real_t rc2 = cutoff() * cutoff();
            logger.debug(std::to_string(cutoff()) + ": Cutoff distance for cell-based particle pairs.");

            // Insert unique particle pairs.
            std::map<std::string, PairList::p_pair_t> pairs{};
            auto insertParticlePair
                = [&pairs] (const p_ptr_t& pi, const p_ptr_t& pj) {
                    // Each particle pair should appear just once in the pair list.
                    auto key = generatePairKey(pi, pj);
                    PairList::p_pair_t ppPair = std::make_pair(pi, pj);
                    auto pair = std::make_pair(key, ppPair);
                    pairs.insert(pair);
                };

            auto neighboring = grid.neighbors();
            for (auto& c : neighboring) {
                // Get central cell.
                auto location = c.first;
                auto central = grid.findCell(location);
                auto centralParticles = central->particles();

                if ( !centralParticles.empty() ) {

                    // Particle pairs in current central cell.
                    for (auto i = 0; i < (centralParticles.size() - 1); ++i) {
                        auto &pi = centralParticles[i];
                        auto ri = pi->position();
                        for (auto j = i + 1; j < centralParticles.size(); ++j) {
                            auto &pj = centralParticles[j];
                            auto rj = pj->position();
                            auto rij = bc->apply(ri, rj);
                            auto rij2 = norm_square<real_t>(rij);
                            if (rij2 <= rc2) {
                                insertParticlePair(pi, pj);
                            }
                        }
                    }
                    auto numberOfCentral = pairs.size();
                    logger.debug(std::to_string(numberOfCentral) +
                                 ": Number of pairs in current central cells.");

                    // Particle pairs between particles in current central cell and its
                    // neighboring cells.
                    auto neighbors = c.second;
                    for (const auto &pi: centralParticles) {
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
                    auto numberOfBetween = pairs.size() - numberOfCentral;
                    logger.debug(std::to_string(numberOfBetween) +
                                 ": Number of pairs between current central and neighboring cells.");
                }

                // Done.
            }
            logger.debug(std::to_string(pairs.size()) +
                         ": Initial number of particle pairs in non-bonded pair list.");

            // Remove bonded particle pairs.
            int counter = 0;
            for (const auto& g: groups) {
                auto& bonds = g->bonds();
                for (const auto& b: bonds) {
                    auto pi = b.getParticleOne();
                    auto pj = b.getParticleTwo();
                    auto key = generatePairKey(pi, pj);
                    pairs.erase(key);
                    counter += 1;
                }
            }
            logger.debug(std::to_string(counter) +
                         ": Number of bonded pairs removed from non-bonded pair list.");

            // Generate final, exclude frozen particle pairs if requested.
            std::size_t nFrozenPairs = 0;
            std::vector<PairList::p_pair_t> particlePairs{};
            for (auto& p: pairs) {
                auto pair = p.second;
                bool frozen = pair.first->frozen() && pair.second->frozen();
                if (!(excludeFrozen && frozen)) {
                    particlePairs.emplace_back(pair);
                } else {
                    nFrozenPairs += 1;
                }
            }
            logger.debug(std::to_string(nFrozenPairs) +
                         ": Number of frozen particle pairs removed from non-bonded pair list.");

            logger.debug(std::to_string(particlePairs.size()) +
            ": Final number of non-bonded particle pairs.");
            auto total = particles.size() * (particles.size() - 1) / 2;
            logger.debug(std::to_string(total) + ": Total number of non-bonded particle pairs.");
            if (total > 0) {
                auto percentage = real_t(particlePairs.size()) * 100.0 / real_t(total);
                logger.debug(std::to_string(percentage) +
                             ": Fraction (%) of total number of non-bonded particle pairs in pair list.");
            } else {
                logger.debug("0: Fraction (%) of total number of non-bonded particle pairs: 0");
            }

            // Done.
            logger.trace("Leaving.");
            return std::move(particlePairs);
        }

        static std::vector<PairList::p_pair_t>
        makePairListConcurrent(const box_ptr_t &box,
                               const dist_t & cutoff,
                               const bc_ptr_t &bc,
                               const std::vector<p_ptr_t>& particles,
                               const std::vector<pg_ptr_t>& groups,
                               const Grid& grid,
                               bool excludeFrozen) {
            static util::Logger logger("simploce::pairlist::makePairListConcurrent()");
            logger.trace("Entering.");

            auto particlePairs =
                    std::move(pairlist::makePairListNonConcurrent(cutoff,
                                                                  bc,
                                                                  particles,
                                                                  groups,
                                                                  grid,
                                                                  excludeFrozen));

            logger.trace("Leaving.");
            return std::move(particlePairs);
        }

    }


    CellPairListGenerator::CellPairListGenerator(param_ptr_t param, bc_ptr_t bc) :
            pair_list_generator{}, param_{std::move(param)}, bc_{std::move(bc)}{
    }

    std::vector<PairList::p_pair_t>
    CellPairListGenerator::generate(const p_system_ptr_t& particleSystem) const {
         static util::Logger logger("simploce::CellPairListGenerator::generate()");
         logger.trace(("Entering."));

         // Setup
         static dist_t cutoff = util::computePairListCutoff(this->param_, particleSystem);
         static Grid grid(particleSystem->box(), bc_, cutoff);
         static bool excludeFrozen = param_->get<bool>("simulation.forces.exclude-frozen");

         auto particlePairs =
             particleSystem->doWithAllFreeGroups<std::vector<PairList::p_pair_t>>([this, particleSystem] (
                 const std::vector<p_ptr_t>& all,
                 const std::vector<p_ptr_t>& free,
                 const std::vector<pg_ptr_t>& groups) {
             grid.assignParticlesToCells(bc_, all);
             return std::move(pairlist::makePairListConcurrent(particleSystem->box(),
                                                               cutoff,
                                                               bc_,
                                                               all,
                                                               groups,
                                                               grid,
                                                               excludeFrozen));
         });

         logger.trace("Leaving.");
         return std::move(particlePairs);
     }

}