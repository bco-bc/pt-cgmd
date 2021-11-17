/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 6, 2019, 12:49 PM
 */

#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/force-field.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/util.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    Interactor::Interactor(sim_param_ptr_t simulationParameters,
                           ff_ptr_t forceField,
                           pair_list_gen_ptr_t pairListGenerator,
                           forces_ptr_t forces,
                           box_ptr_t box,
                           bc_ptr_t bc) :
        simulationParameters_{std::move(simulationParameters)}, forceField_{std::move(forceField)},
        pairListsGenerator_{std::move(pairListGenerator)}, forces_{std::move(forces)},
        box_{std::move(box)}, bc_{std::move(bc)},
        pairLists_{} {
    }

    std::pair<energy_t, energy_t>
    Interactor::interact(const p_system_ptr_t &particleSystem) {
        using energies_t = std::pair<energy_t, energy_t>;

        static util::Logger logger("simploce::Interactor::interact");
        static bool setup = false;
        static std::size_t nUpdatePairLists = 0;
        static std::size_t counter = 0;

        if ( !setup ) {
            nUpdatePairLists = simulationParameters_->get<std::size_t>("simulation.npairlists");
            logger.debug("Number of steps between update of particle pair list: " +
                          util::toString(nUpdatePairLists));
            setup = true;
        }

        // Update particle pair list, if needed.
        if (counter % nUpdatePairLists == 0 || counter == 0) {
            logger.trace("Updating particle pair lists...");
            pairLists_ = std::move(particleSystem->doWithAllFreeGroups<PairLists>([this] (
                    const std::vector<p_ptr_t> &all,
                    const std::vector<p_ptr_t> &free,
                    const std::vector<pg_ptr_t> &groups) {
                return this->pairListsGenerator_->generate(all, free, groups);
            }));

            pairLists_.modified_(true);
            pairLists_.numberOfParticles(particleSystem->numberOfParticles());
            logger.trace("Done.");
        } else {
            pairLists_.modified_(false);
        }

        // Calculate all forces and associated energies.
        particleSystem->resetForces();
        auto energies =
                particleSystem->doWithAllFreeGroups<energies_t>([this] (
                        const std::vector<p_ptr_t>& all,
                        const std::vector<p_ptr_t>& free,
                        const std::vector<pg_ptr_t>& groups) {
                    auto bonded = this->forces_->bonded(all, groups);
                    auto nonBonded = this->forces_->nonBonded(all, this->pairLists_);
                    auto pair = std::make_pair(nonBonded, bonded);
                    return pair;
                });

        counter += 1;

        // Done.
        return energies;
    }

}