/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 6, 2019, 12:49 PM
 */

#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/force-field.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/simulation/forces.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/util.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    Interactor::Interactor(sim_param_ptr_t simulationParameters,
                           pair_list_gen_ptr_t pairListGenerator,
                           forces_ptr_t forces) :
        simulationParameters_{std::move(simulationParameters)},
        pairListsGenerator_{std::move(pairListGenerator)}, forces_{std::move(forces)},
        pairLists_{} {
    }

    std::pair<energy_t, energy_t>
    Interactor::interact(const p_system_ptr_t &particleSystem) {
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
            pairLists_ = pairListsGenerator_->generate(particleSystem);
            pairLists_.modified_(true);
            pairLists_.numberOfParticles(particleSystem->numberOfParticles());
            logger.trace("Done.");
        } else {
            pairLists_.modified_(false);
        }

        // Calculate all forces and associated energies.
        particleSystem->resetForces();
        auto bonded = this->forces_->bonded(particleSystem);
        auto nonBonded = this->forces_->nonBonded(particleSystem, pairLists_);

        // Update counter.
        counter += 1;

        // Done.
        return std::move(std::make_pair(nonBonded, bonded));
    }

    std::pair<energy_t, energy_t>
    Interactor::interact(const p_ptr_t& particle,
                         const p_system_ptr_t &particleSystem) {
        auto energy = this->forces_->interaction(particle, particleSystem);
        return std::make_pair(energy, 0.0);
    }

}