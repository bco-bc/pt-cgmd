/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#include "simploce/potentials/forces.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/potentials/pair-potential.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/potentials/lj.hpp"
#include "simploce/potentials/hp.hpp"
#include "simploce/potentials/halve-attractive-qp.hpp"
#include "simploce/potentials/lj-rf.hpp"
#include "simploce/potentials/rf.hpp"
#include "simploce/potentials/sf.hpp"
#include "simploce/potentials/sc.hpp"
#include "simploce/potentials/hs-sf.hpp"
#include "simploce/potentials/lj-sf.hpp"
#include "simploce/potentials/hs-sc.hpp"
#include "simploce/potentials/wall.hpp"
#include "simploce/potentials/no-pair-potential.hpp"
#include "simploce/potentials/external-potential.hpp"
#include "simploce/potentials/elec-pot-difference.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/bond.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/direction.hpp"
#include "simploce/util/flat-surface.hpp"
#include <utility>
#include <memory>
#include <map>

namespace simploce {

    using result_t = std::pair<energy_t, std::vector<force_t>>;
    using pp_map_t = std::map<std::string, pair_potential_ptr_t>;
    using ep_cont_t = std::vector<std::shared_ptr<external_potential>>;

    static pp_map_t nonBondedPairPotentials_{};
    static pp_map_t bondedPairPotentials_{};

    static ep_cont_t externalPotentials_{};

    static pair_potential_ptr_t
    findPairPotential_(const std::string& key, const pp_map_t& pairPotential) {
        static util::Logger logger("simploce::findPairPotential_()");
        static pair_potential_ptr_t NONE = std::make_shared<NoPairPotential>();
        auto iter = pairPotential.find(key);
        if (iter == pairPotential.end() ) {
            logger.warn(key + ": No pair potential identified.");
            return NONE;
        } else {
            return iter->second;
        }
    }

    /**
     * Associates pair potential with pair of particle specifications.
     * @tparam P Particle typeName.
     * @param forceField Force field.
     * @param box Simulation box.
     * @param bc Boundary condition.
     * @return Associated pair potentials.
     */
    static void
    associatePairPotentials_(const dist_t& cutoff,
                             const std::vector<p_ptr_t>& all,
                             const ff_ptr_t &forceField,
                             const box_ptr_t &box,
                             const bc_ptr_t &bc) {
        static util::Logger logger("simploce::associatePairPotentials_()");

        // Non-bonded pair potentials.
        nonBondedPairPotentials_.clear();
        auto interactionParameters = forceField->nonBondedSpecifications();
        for (auto& ip : interactionParameters) {
            std::string typeName = ip.typeName;
            if (typeName == conf::LJ) {
                std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                auto lj = std::make_shared<LJ>(forceField, bc);
                auto pair1 = std::make_pair(key1, lj);
                nonBondedPairPotentials_.emplace(pair1);
                if ( ip.spec1 != ip.spec2 ) {
                    std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                    auto pair2 = std::make_pair(key2, lj);
                    nonBondedPairPotentials_.emplace(pair2);
                }
            } else if (typeName == conf::LJ_RF) {
                std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                real_t kappa = properties::kappa(all);
                auto rf = std::make_shared<RF>(cutoff, kappa, forceField, box, bc);
                auto lj_rf = std::make_shared<LJ_RF>(kappa, forceField, box, bc, rf);
                auto pair1 = std::make_pair(key1, lj_rf);
                nonBondedPairPotentials_.emplace(pair1);
                if ( ip.spec1 != ip.spec2 ) {
                    std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                    auto pair2 = std::make_pair(key2, lj_rf);
                    nonBondedPairPotentials_.emplace(pair2);
                }
            } else if (typeName == conf::RF) {
                std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                real_t kappa = properties::kappa(all);
                auto rf = std::make_shared<RF>(cutoff, kappa, forceField, box, bc);
                auto pair1 = std::make_pair(key1, rf);
                nonBondedPairPotentials_.emplace(pair1);
                if ( ip.spec1 != ip.spec2 ) {
                    std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                    auto pair2 = std::make_pair(key2, rf);
                    nonBondedPairPotentials_.emplace(pair2);
                }
            } else if (typeName == conf::LJ_SF) {
                std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                auto sf = std::make_shared<SF>(cutoff, forceField, box, bc);
                auto lj_sf = std::make_shared<LJ_SF>(forceField, bc, sf);
                auto pair1 = std::make_pair(key1, lj_sf);
                nonBondedPairPotentials_.emplace(pair1);
                if ( ip.spec1 != ip.spec2 ) {
                    std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                    auto pair2 = std::make_pair(key2, lj_sf);
                    nonBondedPairPotentials_.emplace(pair2);
                }
            } else if (typeName == conf::HS_SF ) {
               std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
               auto sf = std::make_shared<SF>(cutoff, forceField, box, bc);
               auto hs_sf = std::make_shared<HS_SF>(forceField, bc, sf);
               auto pair1 = std::make_pair(key1, hs_sf);
               nonBondedPairPotentials_.emplace(pair1);
               if ( ip.spec1 != ip.spec2 ) {
                    std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                    auto pair2 = std::make_pair(key2, hs_sf);
                    nonBondedPairPotentials_.emplace(pair2);
                }
            } else if (typeName == conf::HS_SC) {
                std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                auto hs_sc = std::make_shared<HS_SC>(forceField, bc);
                auto pair1 = std::make_pair(key1, hs_sc);
                nonBondedPairPotentials_.emplace(pair1);
                if ( ip.spec1 != ip.spec2 ) {
                    std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                    auto pair2 = std::make_pair(key2, hs_sc);
                    nonBondedPairPotentials_.emplace(pair2);
                }
            } else {
                util::logAndThrow(logger, "Pair potential '" + typeName + "' is not available.");
            }
        }

        // Bonded pair potentials.
        bondedPairPotentials_.clear();
        interactionParameters = forceField->bondedSpecifications();
        for (auto& ip : interactionParameters) {
            std::string type = ip.typeName;
            if (type == conf::HP ) {
                std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                auto hp = std::make_shared<HP>(forceField, bc);
                auto pair1 = std::make_pair(key1, hp);
                bondedPairPotentials_.emplace(pair1);
                if ( ip.spec1 != ip.spec2 ) {
                    std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                    auto pair2 = std::make_pair(key2, hp);
                    bondedPairPotentials_.emplace(pair2);
                }
            } else if ( type == conf::HA_QP ) {
                std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                auto ha_qp = std::make_shared<HalveAttractiveQP>(forceField, bc);
                auto pair1 = std::make_pair(key1, ha_qp);
                bondedPairPotentials_.emplace(pair1);
                if (ip.spec1 != ip.spec2) {
                    std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                    auto pair2 = std::make_pair(key2, ha_qp);
                    bondedPairPotentials_.emplace(pair2);
                }
            } else {
                util::logAndThrow(logger,"Pair potential '" + type + "' is not available.");
            }
        }

        // Log some info.
        logger.debug("Number of associated non-bonded pair potentials: " +
                      util::toString(nonBondedPairPotentials_.size()));
        logger.debug("Number of associated bonded pair potentials: " +
                      util::toString(bondedPairPotentials_.size()));
    }

    static void
    associateExternalPotentials_(const ff_ptr_t &forceField,
                                 const bc_ptr_t &bc) {
        util::Logger logger("simploce::associateExternalPotentials_()");
        auto external = forceField->externalSpecifications();
        for (const auto& spec : external) {
            std::string typeName = spec.typeName;
            if (typeName == conf::ELECTRIC_POTENTIAL_DIFFERENCE) {
                auto direction = Direction::valueOf(spec.direction);
                auto potential = std::make_shared<ElectricPotentialDifference>(spec.deltaV,
                                                                               spec.distance,
                                                                               spec.eps_inside_rc,
                                                                               bc,
                                                                               direction);
                externalPotentials_.emplace_back(potential);
            } else if (typeName == conf::WALL) {
                auto plane = Plane::valueOf(spec.plane);
                dist_t distance{spec.distance};
                FlatSurface flatSurface{plane, distance};
                srf_charge_density_t sigma{spec.sigma};
                auto potential = std::make_shared<Wall>(spec.C12, spec.C6, bc, flatSurface, sigma);
                externalPotentials_.emplace_back(potential);
            } else {
                util::logAndThrow(logger,"External potential '" + typeName + "' is not available.");
            }
        }

        // Log some info.
        logger.debug("Number of external potentials: " + util::toString(externalPotentials_.size()));
    }

    /**
     * Returns non-bonded potential energy forces for given pairs of particles due to interactions
     * between these particles.
     * @return Potential energy and forces on particles.
     */
    static std::pair<energy_t, std::vector<force_t>>
    computeNonBondedForcesParticlePairs_(const PairLists::pp_pair_cont_t &particlePairs,
                                         const pp_map_t &nonBondedPairPotentials,
                                         int numberOfParticles) {
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
           pair_potential &pairPotential = *findPairPotential_(key, nonBondedPairPotentials);
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
     * Computes non-bonded forces and associated interaction potential energies.
     * @return Non-bonded potential energy.
     */
    static energy_t
    computeNonBondedForces_(const std::vector<p_ptr_t> &all,
                            const PairLists &pairLists,
                            const pp_map_t &nonBondedPairPotentials,
                            const box_ptr_t &box,
                            const bc_ptr_t &bc) {
        using pp_pair_cont_t = PairLists::pp_pair_cont_t;

        static util::Logger logger("simploce::forces::computeNonBondedForces_()");
        static std::vector<pp_pair_cont_t> subPairLists;
        static bool firstTime = true;

        if ( pairLists.isEmpty() ) {
            return 0.0;
        }

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
                                        computeNonBondedForcesParticlePairs_,
                                        std::ref(single),
                                        std::ref(nonBondedPairPotentials),
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
                    auto result = computeNonBondedForcesParticlePairs_(single, nonBondedPairPotentials,
                                                                       int(numberOfParticles));
                    results.push_back(result);
                }
            } else {
                // Sequentially. All particle-particle pair interactions are handled by the
                // current thread.
                if ( firstTime ) {
                    logger.debug("Lennard-Jones and electrostatic forces and energies are computed sequentially.");
                }
                auto result = computeNonBondedForcesParticlePairs_(pairLists.particlePairList(),
                                                                   nonBondedPairPotentials,
                                                                   int(numberOfParticles));
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
     * Computes bonded forces and associated potential energies. Forces are added to particles.
     * @return Bonded interaction potential energy.
     */
    static energy_t
    computeBondedForces_(const std::vector<pg_ptr_t> &groups,
                         const pp_map_t &bondedPairPotentials) {
        energy_t bonded{0.0};
        for (auto &group : groups) {
            for (auto &bond : group->bonds() ) {
                // Get bonded particles.
                auto p1 = bond.getParticleOne();
                auto p2 = bond.getParticleTwo();

                // Call pair potential.
                std::string key = p1->spec()->name() + "-" + p2->spec()->name();
                pair_potential& pairPotential = *findPairPotential_(key, bondedPairPotentials);
                std::pair<energy_t, force_t> ef = pairPotential(p1, p2);

                // Store results.
                bonded += std::get<0>(ef);
                auto f = std::get<1>(ef);
                auto f1 = p1->force();
                f1 += f;
                p1->force(f1);
                auto f2 = p2->force();
                f2 -= f;
                p2->force(f2);
            }
        }
        return bonded;
    }

    static energy_t
    computeExternalForces(const std::vector<p_ptr_t>& all,
                          const ep_cont_t& externalPotentials) {
        energy_t external{};
        for (auto& particle: all) {
            for (auto& ep: externalPotentials) {
                auto result = (*ep)(particle);
                external += result.first;
                auto f = particle->force();
                f += result.second;
                particle->force(f);
            }
        }
        return external;
    }

    /**
     * Returns interaction energy of a given particle with all other particles (except with itself)
     * within the cutoff distance.
     * @param particle Particle.
     * @param all -All- particles.
     * @param nonBondedPairPotentials Pair potentials.
     * @return Bonded, non-bonded, and external potential energy.
     */
    static std::tuple<energy_t, energy_t, energy_t>
    interactionEnergyOfOneParticle_(const p_ptr_t& particle,
                                    const std::vector<p_ptr_t>& free,
                                    const std::vector<pg_ptr_t>& groups,
                                    const pp_map_t& nonBondedPairPotentials,
                                    const pp_map_t& bondedPairPotentials,
                                    const ep_cont_t& externalPotentials,
                                    const box_ptr_t& box,
                                    const dist_t cutoff,
                                    const bc_ptr_t& bc) {
        static auto rc2 = cutoff() * cutoff();

        energy_t nonBonded{0.0}, bonded{0.0}, external{0.0};
        auto name = particle->spec()->name();
        auto id = particle->id();
        auto pr = particle->position();

        // Non-bonded interaction with free particles.
        for (const auto& p: free) {
            if ( p->id() != id) {
                auto r = p->position();
                dist_vect_t R = bc->apply(pr, r);
                auto R2 = norm_square<real_t>(R);
                if ( R2 <= rc2 ) {
                    std::string key = name + "-" + p->spec()->name();
                    auto &pairPotential = *nonBondedPairPotentials.at(key);
                    // Force is ignored.
                    std::pair<energy_t, force_t> ef = pairPotential(particle, p);
                    nonBonded += std::get<0>(ef);
                }
            }
        }

        // Interaction with and within groups.
        for (const auto& g: groups) {
            if (!g->contains(particle)) {
                auto particles = g->particles();
                for (const auto& p: particles) {
                    auto r = p->position();
                    dist_vect_t R = bc->apply(pr, r);
                    auto R2 = norm_square<real_t>(R);
                    if ( R2 <= rc2 ) {
                        std::string key = name + "-" + p->spec()->name();
                        auto &pairPotential = *nonBondedPairPotentials.at(key);
                        // Force is ignored.
                        std::pair<energy_t, force_t> ef = pairPotential(particle, p);
                        nonBonded += std::get<0>(ef);
                    }
                }
            } else {
                auto particles = g->particles();
                for (const auto& p: particles) {
                    if (particle->id() != p->id() ) {
                        auto r = p->position();
                        dist_vect_t R = bc->apply(pr, r);
                        std::string key = name + "-" + p->spec()->name();
                        auto &pairPotential = *bondedPairPotentials.at(key);
                        // Force is ignored.
                        std::pair<energy_t, force_t> ef = pairPotential(particle, p);
                        bonded += std::get<0>(ef);
                    }
                }
            }
        }

        // Interaction with external potentials.
        for (auto& ep: externalPotentials) {
            auto& potential = *ep;
            std::pair<energy_t, force_t> ef = potential(particle);
            // Force is ignored.
            external += std::get<0>(ef);
        }
        return std::move(std::make_tuple(bonded, nonBonded, external));
    }

    static void
    associatePotentials_(const dist_t & cutoff,
                         const std::vector<p_ptr_t>& all,
                         const ff_ptr_t& forceField,
                         const box_ptr_t& box,
                         const bc_ptr_t& bc) {
        static bool ASSOCIATED{false};
        if (ASSOCIATED) {
            return;
        }
        associatePairPotentials_(cutoff, all, forceField, box, bc);
        associateExternalPotentials_(forceField, bc);
        ASSOCIATED = true;
    }

    Forces::Forces(const param_ptr_t& param, bc_ptr_t bc, ff_ptr_t forceField) :
            cutoff_{}, bc_{std::move(bc)}, forceField_{std::move(forceField)} {
        cutoff_ = param->get<real_t>("forces.nb.cutoff");
    }

    energy_t Forces::nonBonded(const p_system_ptr_t& particleSystem,
                               const PairLists &pairLists) {
        auto box = particleSystem->box();
        return particleSystem->doWithAll<energy_t>([this, pairLists, box] (std::vector<p_ptr_t>& all) {
            associatePotentials_(cutoff_, all, forceField_, box, this->bc_);
            return computeNonBondedForces_(all,
                                           pairLists,
                                           nonBondedPairPotentials_,
                                           box,
                                           this->bc_);
        });
    }

    energy_t
    Forces::external(const p_system_ptr_t& particleSystem) {
        auto box = particleSystem->box();
        return particleSystem->doWithAll<energy_t>([this, box] (const std::vector<p_ptr_t>& all) {
            associatePotentials_(cutoff_, all, this->forceField_, box, this->bc_);
            return computeExternalForces(all, externalPotentials_);
        });
    }

    energy_t
    Forces::bonded(const p_system_ptr_t& particleSystem) {
        auto box = particleSystem->box();
        return particleSystem->doWithAllFreeGroups<energy_t>([this, box] (
                const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups
                ) {
            associatePotentials_(cutoff_, all, this->forceField_, box, this->bc_);
            return computeBondedForces_(groups, bondedPairPotentials_);
        });
    }

    std::tuple<energy_t, energy_t, energy_t>
    Forces::interaction(const p_ptr_t& particle,
                        const p_system_ptr_t& particleSystem) {
        auto box = particleSystem->box();
        return particleSystem->doWithAllFreeGroups<std::tuple<energy_t, energy_t, energy_t>>([this, particle, box] (
                const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) {
            associatePotentials_(cutoff_, all, this->forceField_, box, this->bc_);
            return interactionEnergyOfOneParticle_(particle,
                                                   free,
                                                   groups,
                                                   nonBondedPairPotentials_,
                                                   bondedPairPotentials_,
                                                   externalPotentials_,
                                                   box,
                                                   this->cutoff_,
                                                   this->bc_);
        });
    }

}
