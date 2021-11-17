/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#include "simploce/simulation/forces.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/pair-potential.hpp"
#include "simploce/simulation/lj.hpp"
#include "simploce/simulation/hp.hpp"
#include "simploce/simulation/halve-attractive-qp.hpp"
#include "simploce/simulation/lj-rf.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/util/util.hpp"
#include <utility>
#include <memory>
#include <map>

namespace simploce {

     using result_t = std::pair<energy_t, std::vector<force_t>>;

     static std::map<std::string, std::shared_ptr<pair_potential<Atom>>> aaPairPotentials{};

     static std::map<std::string, std::shared_ptr<pair_potential<Bead>>> bbPairPotentials{};

     /**
      * Associates pair potential with pair of particle specifications.
      * @tparam P Particle type.
      * @param forceField Force field.
      * @param box Simulation box.
      * @param bc Boundary condition.
      * @return Associated pair potentials.
      */
     template <typename P>
     static std::map<std::string, std::shared_ptr<pair_potential<P>>>
     associatePairPotentials_(const std::vector<std::shared_ptr<P>>& all,
                              const ff_ptr_t &forceField,
                              const box_ptr_t &box,
                              const bc_ptr_t &bc) {
         using pair_pot_ptr_t = std::shared_ptr<pair_potential<P>>;

         static util::Logger logger("simploce::associatePairPotentials_()");

         std::map<std::string, pair_pot_ptr_t> pairPotentials{};
         auto& interactionParameters = forceField->interactionSpecifications();
         for (auto& ip : interactionParameters) {
             std::string type = ip.type;
             if (type == conf::LJ) {
                 std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                 auto lj = std::make_shared<LJ<P>>(forceField, box, bc);
                 auto pair1 = std::make_pair(key1, lj);
                 pairPotentials.emplace(pair1);
                 if ( ip.spec1 != ip.spec2 ) {
                     std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                     auto pair2 = std::make_pair(key2, lj);
                     pairPotentials.emplace(pair2);
                 }
             } else if (type == conf::LJ_RF) {
                 std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                 real_t kappa = properties::kappa(all);
                 auto lj_rf = std::make_shared<LJ_RF<P>>(kappa, forceField, box, bc);
                 auto pair1 = std::make_pair(key1, lj_rf);
                 pairPotentials.emplace(pair1);
                 if ( ip.spec1 != ip.spec2 ) {
                     std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                     auto pair2 = std::make_pair(key2, lj_rf);
                     pairPotentials.emplace(pair2);
                 }
             }
             else if (type == conf::HP ) {
                 std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                 auto hp = std::make_shared<HP<P>>(forceField, bc);
                 auto pair1 = std::make_pair(key1, hp);
                 pairPotentials.emplace(pair1);
                 if ( ip.spec1 != ip.spec2 ) {
                     std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                     auto pair2 = std::make_pair(key2, hp);
                     pairPotentials.emplace(pair2);
                 }
             } else if ( type == conf::HA_QP ) {
                 std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                 auto ha_qp = std::make_shared<HalveAttractiveQP<P>>(forceField, bc);
                 auto pair1 = std::make_pair(key1, ha_qp);
                 pairPotentials.emplace(pair1);
                 if ( ip.spec1 != ip.spec2 ) {
                     std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                     auto pair2 = std::make_pair(key2, ha_qp);
                     pairPotentials.emplace(pair2);
                 }
             }
         }

         // Log some info.
         logger.debug("Number of assigned pair potentials: " + util::toString(pairPotentials.size()));

         return std::move(pairPotentials);
     }

    /**
     * Returns (non-bonded) interaction forces for given pairs of particles.
     */
    template <typename P>
    static std::pair<energy_t, std::vector<force_t>>
    computeParticlePairForces_(const typename PairLists<P>::pp_pair_cont_t &particlePairs,
                               const std::map<std::string, std::shared_ptr<pair_potential<P>>> &pairPotentials,
                               int numberOfParticles) {
        // Particle pointer type.
        using p_ptr_t = std::shared_ptr<P>;

        // Compute interactions for all particle pairs.
        std::vector<force_t> forces(numberOfParticles, force_t{0.0, 0.0, 0.0});
        energy_t energy{0.0};

        if ( particlePairs.empty() ) {
            return std::move(std::make_pair(energy, forces));
        }

        for (auto &particlePair : particlePairs) {

            // First particle
            p_ptr_t pi = particlePair.first;
            auto index_i = pi->index();

            // Second particle
            p_ptr_t pj = particlePair.second;
            auto index_j = pj->index();

            // Call pair potential.
            std::string key = pi->spec()->name() + "-" + pj->spec()->name();
            pair_potential<P>& pairPotential = *pairPotentials.at(key);
            std::pair<energy_t, force_t> ef = pairPotential(pi, pj);

            // Store results.
            energy += std::get<0>(ef);
            forces[index_i] += std::get<1>(ef);
            forces[index_j] -= std::get<1>(ef);
        }

        // Done.
        return std::make_pair(energy, forces);
    }

    /**
     * Computes non-bonded forces and associated potential energies.
     */
    template <typename P>
    static energy_t
    computeNonBondedForces_(const std::vector<std::shared_ptr<P>> &all,
                            const PairLists<P> &pairLists,
                            const std::map<std::string, std::shared_ptr<pair_potential<P>>> &pairPotentials,
                            const box_ptr_t &box,
                            const bc_ptr_t &bc) {

            using pp_pair_cont_t = typename PairLists<P>::pp_pair_cont_t;

            if ( pairLists.isEmpty() ) {
                return 0.0;
            }

            static util::Logger logger("simploce::forces::computeNonBondedForces_()");
            static std::vector<pp_pair_cont_t> subPairLists;
            static bool firstTime = true;

            auto numberOfParticles = pairLists.numberOfParticles();
            std::vector<result_t> results{};

            if ( numberOfParticles > simploce::conf::MIN_NUMBER_OF_PARTICLES ) {
                // Handle particle-particle interaction concurrently as tasks, where one task is executed
                // by the current thread.
                if ( firstTime ) {
                    logger.debug("Lennard-Jones and electrostatic forces and energies are computed concurrently.");
                }
                std::vector<std::future<result_t> > futures{};
                if ( pairLists.isModified() || firstTime) {
                    subPairLists = util::makeSubLists(pairLists.particlePairList());
                }
                auto numberOfTasks = subPairLists.size() - 1;
                if ( numberOfTasks > 1 ) {
                    // Set up concurrent force calculations.
                    for (std::size_t k = 0; k != numberOfTasks; ++k) {
                        const pp_pair_cont_t& single = subPairLists[k];
                        futures.push_back(
                                std::async(
                                        std::launch::async,
                                        computeParticlePairForces_<P>,
                                        std::ref(single),
                                        std::ref(pairPotentials),
                                        numberOfParticles
                                )
                        );
                    }
                    // Wait for tasks to complete.
                    results = util::waitForAll<result_t>(futures);
                }
                // One remaining set of particle-particle interactions is handled
                // by the current thread.
                const pp_pair_cont_t &single = *(subPairLists.end() - 1);
                if ( !single.empty() ) {
                    auto result = computeParticlePairForces_<P>(single, pairPotentials, numberOfParticles);
                    results.push_back(result);
                }
            } else {
                // Sequentially. All particle-particle pair interactions are handled by the
                // current thread.
                if ( firstTime ) {
                    logger.debug("Lennard-Jones and electrostatic forces and energies are computed sequentially.");
                }
                auto result = computeParticlePairForces_<P>(pairLists.particlePairList(), pairPotentials,
                                                            numberOfParticles);
                results.push_back(result);
            }

            // Collect potential energies and forces.
            energy_t energy{0.0};
            for (auto& result : results) {
                const auto& forces = result.second;
                for (auto& particle : all) {
                    auto index = particle->index();
                    force_t f = forces[index] + particle->force();
                    particle->force(f);
                }
                energy += result.first;
            }

            // Done. Returns potential energy.
            firstTime = false;
            return energy;
    }


    /**
     * Computes bonded forces and associated potential energies.
     */
    template <typename P>
    static energy_t
    computeBondedForces(const std::vector<std::shared_ptr<ParticleGroup<P>>> &groups,
                        const std::map<std::string, std::shared_ptr<pair_potential<P>>> &pairPotentials) {
        energy_t energy{0.0};
        for (auto &group : groups) {
            for (auto &bond : group->bonds() ) {
                // Get bonded particles.
                auto& p1 = bond.getParticleOne();
                auto index_1 = p1->index();
                auto& p2 = bond.getParticleTwo();
                auto index_2 = p2->index();

                // Call pair potential.
                std::string key = p1->spec()->name() + "-" + p2->spec()->name();
                pair_potential<P>& pairPotential = *pairPotentials.at(key);
                std::pair<energy_t, force_t> ef = pairPotential(p1, p2);

                // Store results.
                energy += std::get<0>(ef);
                auto f = std::get<1>(ef);
                auto f1 = p1->force();
                f1 += f;
                p1->force(f1);
                auto f2 = p2->force();
                f2 -= f;
                p2->force(f2);
            }
        }
        return energy;
    }

    Forces::Forces(box_ptr_t box, bc_ptr_t bc, ff_ptr_t forceField) :
        box_{std::move(box)}, bc_{std::move(bc)}, forceField_{std::move(forceField)} {
    }

    energy_t Forces::nonBonded(const std::vector<atom_ptr_t> &all,
                               const PairLists<Atom> &pairLists) {
        if ( aaPairPotentials.empty() ) {
            aaPairPotentials = associatePairPotentials_<Atom>(all, forceField_, box_, bc_);
        }
        return computeNonBondedForces_<Atom>(all, pairLists, aaPairPotentials, box_, bc_);
    }

    energy_t Forces::bonded(const std::vector<atom_ptr_t> &all,
                            const std::vector<std::shared_ptr<ParticleGroup<Atom>>> &groups) {
        if ( aaPairPotentials.empty() ) {
            aaPairPotentials = associatePairPotentials_<Atom>(all, forceField_, box_, bc_);
        }
        return computeBondedForces<Atom>(groups, aaPairPotentials);
    }

    energy_t Forces::nonBonded(const std::vector<bead_ptr_t> &all,
                               const PairLists<Bead> &pairLists) {
        if ( bbPairPotentials.empty() ) {
            bbPairPotentials = associatePairPotentials_<Bead>(all, forceField_, box_, bc_);
        }
        return computeNonBondedForces_<Bead>(all, pairLists, bbPairPotentials, box_, bc_);
    }

    energy_t Forces::bonded(const std::vector<bead_ptr_t> &all,
                            const std::vector<std::shared_ptr<ParticleGroup<Bead>>> &groups) {
        if ( bbPairPotentials.empty() ) {
            bbPairPotentials = associatePairPotentials_<Bead>(all, forceField_, box_, bc_);
        }
        return computeBondedForces<Bead>(groups, bbPairPotentials);
    }
}