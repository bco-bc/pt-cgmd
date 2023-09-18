/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 6, 2019, 12:49 PM
 */

#include "simploce/simulation/interactor.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/potentials/forces.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/bond.hpp"
#include "simploce/util/util.hpp"
#include <memory>
#include <utility>
#include <vector>

namespace simploce {
    namespace interactor {

        using p_pair_t = PairList::p_pair_t;

        static std::vector<p_pair_t>
        bondedParticlePairs(const p_system_ptr_t &particleSystem) {
            auto particlePairs = particleSystem->doWithAllFreeGroups<std::vector<p_pair_t>>([](
                std::vector<p_ptr_t> &all, std::vector<p_ptr_t> &free, std::vector<pg_ptr_t> &groups) {
                std::vector<p_pair_t> particlePairs{};
                for (const auto &g: groups) {
                    auto bonds = g->bonds();
                    for (auto &b: bonds) {
                        p_pair_t pair = std::make_pair(b.getParticleOne(), b.getParticleTwo());
                        particlePairs.emplace_back(pair);
                    }
                }
                return std::move(particlePairs);
            });
            return std::move(particlePairs);
        }
    }
    
    Interactor::Interactor(param_ptr_t param,
                           pair_list_gen_ptr_t pairListGenerator,
                           forces_ptr_t forces) :
            param_{std::move(param)}, pairListGenerator_{std::move(pairListGenerator)},
            forces_{std::move(forces)}, pairList_{std::make_shared<PairList>()} {
    }

    std::tuple<energy_t, energy_t, energy_t>
    Interactor::interact(const p_system_ptr_t &particleSystem) {
        static util::Logger logger("simploce::Interactor::interact()");
        logger.trace("Entering.");

        static auto nUpdatePairLists = param_->get<std::size_t>("simulation.npairlists");
        static auto includeExternal = param_->get<bool>("simulation.forces.include-external");
        static std::size_t counter = 0;

        // Set the bonded pair list. Only once.
        if (counter == 0) {
            logger.debug(std::to_string(nUpdatePairLists) +
                         ": Number of steps between an update of the particle pair list.");
            auto particlePairs = interactor::bondedParticlePairs(particleSystem);
            pairList_->bondedParticlePairs(particlePairs);
            logger.debug(std::to_string(pairList_->bondedParticlePairs().size()) +
                         ": Number of bonded particle pairs.");
        }

        // Update particle pair list, if needed.
        if (counter % nUpdatePairLists == 0 || counter == 0) {
            auto particlePairs = pairListGenerator_->generate(particleSystem);
            pairList_->nonBoundedParticlePairs(particlePairs);
            pairList_->modified(true);
            pairList_->numberOfParticles(particleSystem->numberOfParticles());
        } else {
            pairList_->modified(false);
        }

        // Calculate all forces and associated energies.
        particleSystem->resetForces();
        auto bonded = forces_->bonded(particleSystem);
        auto nonBonded = forces_->nonBonded(particleSystem, pairList_);
        energy_t external{0};
        if (includeExternal) {
            if (counter == 0)
                logger.warn("Including external potentials.");
            external = this->forces_->external(particleSystem);
        }

        // Update counter.
        counter += 1;

        // Done.
        logger.trace("Leaving.");
        return std::move(std::tuple<energy_t, energy_t, energy_t>{bonded, nonBonded, external});
    }

    std::tuple<energy_t, energy_t, energy_t>
    Interactor::interact(const p_ptr_t& particle,
                         const p_system_ptr_t &particleSystem) {
        return std::move(this->forces_->interaction(particle, particleSystem));
    }

    pairlist_ptr_t
    Interactor::pairList() const {
        return pairList_;
    }

}