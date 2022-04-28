/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 6, 2019, 12:49 PM
 */

#include "simploce/simulation/interactor.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/potentials/forces.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/util.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    Interactor::Interactor(param_ptr_t param,
                           pair_list_gen_ptr_t pairListGenerator,
                           forces_ptr_t forces) :
            param_{std::move(param)},
            pairListsGenerator_{std::move(pairListGenerator)}, forces_{std::move(forces)},
            pairLists_{} {
    }

    std::tuple<energy_t, energy_t, energy_t>
    Interactor::interact(const p_system_ptr_t &particleSystem) {
        static util::Logger logger("simploce::Interactor::interact");

        static auto nUpdatePairLists = param_->get<std::size_t>("simulation.npairlists");
        static auto includeExternal = param_->get<bool>("simulation.include-external");
        static std::size_t counter = 0;

        // Update particle pair list, if needed.
        if (counter % nUpdatePairLists == 0 || counter == 0) {
            logger.debug("Number of steps between update of particle pair list: " +
                          util::toString(nUpdatePairLists));
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
        energy_t external{0};
        if (includeExternal) {
            external = this->forces_->external(particleSystem);
        }

        // Update counter.
        counter += 1;

        // Done.
        return std::move(std::tuple<energy_t, energy_t, energy_t>{bonded, nonBonded, external});
    }

    std::tuple<energy_t, energy_t, energy_t>
    Interactor::interact(const p_ptr_t& particle,
                         const p_system_ptr_t &particleSystem) {
        return std::move(this->forces_->interaction(particle, particleSystem));
    }

}