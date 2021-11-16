/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 12:31 PM
 */

#ifndef INTERACTOR_HPP
#define INTERACTOR_HPP

#include "force-field.hpp"
#include "sim-data.hpp"
#include "pair-lists.hpp"
#include "s-types.hpp"
#include "forces.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/util/logger.hpp"
#include <memory>
#include <utility>
#include <vector>
#include <utility>
#include <string>

namespace simploce {

    /**
     * Creates atom pair lists.
     * @param pairListGenerator Pair lists generator.
     * @param atomistic Atomistic particle system.
     * @return PairLists.
     */
    PairLists<Atom>
    updatePairLists(const std::shared_ptr<pair_lists_generator<Atom>> &pairListGenerator,
                    const std::shared_ptr<ParticleSystem<Atom, ParticleGroup<Atom>>> &atomistic);

    /**
     * Creates bead pair lists.
     * @param pairListGenerator Pair lists generator.
     * @param coarseGrained Coarse grained particle system.
     * @return Pair lists.
     */
    PairLists<Bead>
    updatePairLists(const std::shared_ptr<pair_lists_generator<Bead>> &pairListGenerator,
                    const std::shared_ptr<ParticleSystem<Bead, ParticleGroup<Bead>>> &coarseGrained);
    
    /**
     * "One that interacts". Its responsibility is ensure that the particle pair list
     * is (re)evaluated on a regular basis. Its sits in between a simulation and the
     * force calculation.
     * @param P Particle type.
     */
    template <typename P>
    class Interactor {
    public:

        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;

        /**
         * Particle group point type.
         */
         using pg_ptr_t = std::shared_ptr<ParticleGroup<P>>;

        /**
         * Particle system pointer type.
         */
        using p_sys_ptr_t = std::shared_ptr<ParticleSystem<P, ParticleGroup<P>>>;

        /**
         * Particle pair list generator pointer type.
         */
        using pair_list_gen_ptr_t = std::shared_ptr<pair_lists_generator<P>>;

        /**
         * Constructor. All arguments are required.
         * @param simulationParameters Simulation parameters.
         * @param forceField Force field.
         * @param PairListGenerator Particle pair list generator.
         * @param box Simulation box.
         */
        Interactor(sim_param_ptr_t simulationParameters,
                   ff_ptr_t forceField,
                   pair_list_gen_ptr_t pairListGenerator,
                   box_ptr_t box,
                   bc_ptr_t bc);

        /**
         * Calculates forces on all particles in the given particle system.
         * @param simulationParameters Simulation parameters,
         * @param particleSystem Particle system.
         * @return Returns non-bonded and bonded potential energy, respectively.
         */
        std::pair<energy_t, energy_t> interact(const p_sys_ptr_t &particleSystem);

    private:

        sim_param_ptr_t simulationParameters_;
        ff_ptr_t forceField_;
        pair_list_gen_ptr_t pairListsGenerator_;
        box_ptr_t box_;
        bc_ptr_t bc_;

        PairLists<P> pairLists_;
        Forces forces_;

    };

    template <typename P>
    Interactor<P>::Interactor(sim_param_ptr_t simulationParameters,
                              ff_ptr_t forceField,
                              pair_list_gen_ptr_t pairListGenerator,
                              box_ptr_t box,
                              bc_ptr_t bc) :
        simulationParameters_{std::move(simulationParameters)}, forceField_{std::move(forceField)},
        pairListsGenerator_{std::move(pairListGenerator)}, box_{std::move(box)}, bc_{std::move(bc)},
        pairLists_{}, forces_{box_, bc_, forceField_} {
    }

    template <typename P>
    std::pair<energy_t, energy_t> Interactor<P>::interact(const p_sys_ptr_t &particleSystem) {
        // Return type.
        using energies_t = std::pair<energy_t, energy_t>;

        static util::Logger logger("simploce::Interactor<P>::interact");

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
            pairLists_ = updatePairLists(pairListsGenerator_, particleSystem);
            pairLists_.modified_(true);
            pairLists_.numberOfParticles(particleSystem->numberOfParticles());
            logger.trace("Done.");
        } else {
            pairLists_.modified_(false);
        }

        // Calculate all forces and associated energies.
        particleSystem->resetForces();
        auto energies =
                particleSystem->template doWithAllFreeGroups<energies_t>([this] (
                        const std::vector<p_ptr_t>& all,
                        const std::vector<p_ptr_t>& free,
                        const std::vector<pg_ptr_t>& groups) {
                    auto bonded = this->forces_.bonded(all, groups);
                    logger.debug("Bonded energy: " + util::toString(bonded));
                    auto nonBonded = this->forces_.nonBonded(all, this->pairLists_);
                    logger.debug("Non-bonded energy: " + util::toString(nonBonded));
                    return std::make_pair(nonBonded, bonded);
                });

        counter += 1;

        // Done.
        return energies;
    }

}

#endif /* INTERACTOR_HPP */

