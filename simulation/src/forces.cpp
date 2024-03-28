/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#include "simploce/potentials/forces.hpp"
#include "simploce/simulation/pair-list.hpp"
#include "simploce/potentials/cutoffs.hpp"
#include "simploce/conf/sp-conf.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/potentials/pair-potential.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/potentials/lj.hpp"
#include "simploce/potentials/hp.hpp"
#include "simploce/potentials/halve-attractive-qp.hpp"
#include "simploce/potentials/halve-attractive-hp.hpp"
#include "simploce/potentials/lj-rf.hpp"
#include "simploce/potentials/rf.hpp"
#include "simploce/potentials/sf.hpp"
#include "simploce/potentials/sc.hpp"
#include "simploce/potentials/hs-sf.hpp"
#include "simploce/potentials/lj-sf.hpp"
#include "simploce/potentials/hs-sc.hpp"
#include "simploce/potentials/soft-repulsion.hpp"
#include "simploce/potentials/hp-sr.hpp"
#include "simploce/potentials/gauss-sf.hpp"
#include "simploce/potentials/gauss-sf-sr.hpp"
#include "simploce/potentials/lekner.hpp"
#include "simploce/potentials/hs-lekner.hpp"
#include "simploce/potentials/wall.hpp"
#include "simploce/potentials/non-interacting.hpp"
#include "simploce/potentials/external-potential.hpp"
#include "simploce/potentials/voltage.hpp"
#include "simploce/potentials/pressure-gradient.hpp"
#include "simploce/potentials/uniform-surface-charge-density.hpp"
#include "simploce/potentials/vplanes.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/bond.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/flat-surface.hpp"
#include <utility>
#include <memory>
#include <map>

namespace simploce {
    namespace forces {

        using result_t = std::pair<energy_t, std::vector<force_t>>;
        using pp_map_t = std::map<std::string, pair_potential_ptr_t>;
        using ep_cont_t = std::vector<std::shared_ptr<external_potential>>;

        static pp_map_t nonBondedPairPotentials_{};
        static pp_map_t bondedPairPotentials_{};

        static ep_cont_t externalPotentials_{};

        /**
         * Selects pair potential.
         * @param key Key. Identifies pair potential.
         * @param pairPotentials Available pair potentials.
         * @return Selected pair potential. NonInteracting is returned if key cannot be associated with any pair
         * potential.
         * @throws if key cannot be associated with any pair potential.
         */
        static pair_potential_ptr_t
        selectPairPotential(const std::string &key, const pp_map_t &pairPotentials) {
            static util::Logger logger("simploce::forces::selectPairPotential()");
            auto iter = pairPotentials.find(key);
            if (iter == pairPotentials.end()) {
                std::string msg = key + ": No pair potential identified.";
                util::logAndThrow(logger, msg);
                return nullptr;
            } else {
                return iter->second;
            }
        }

        /**
         * Associates pair potential with pair of particle specifications.
         * @tparam P Particle typeName.
         * @param forceField Force field.
         * @param mesoscopic Forces are for a mesoscopic simulation, e.g. DPD.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Associated pair potentials.
         */
        static void
        associatePairPotentials(const rc_ptr_t &cutoffs,
                                const std::vector<p_ptr_t> &all,
                                const ff_ptr_t &forceField,
                                bool mesoscopic,
                                const box_ptr_t &box,
                                const bc_ptr_t &bc) {
            static util::Logger logger("simploce::forces::associatePairPotentials()");
            logger.trace("Entering.");

            // Non-bonded pair potentials.
            auto cutoffSR = cutoffs->shortRanged();
            auto cutoffLR = cutoffs->longRanged();
            nonBondedPairPotentials_.clear();
            auto interactionParameters = forceField->nonBondedSpecifications();
            for (auto &ip: interactionParameters) {
                std::string typeName = ip.typeName;
                logger.debug(typeName + ": Non-bonded pair potential.");
                if (typeName == conf::LJ) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto lj = std::make_shared<LJ>(forceField, bc);
                    auto pair1 = std::make_pair(key1, lj);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, lj);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::LJ_RF) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    real_t kappa = properties::kappa(all);
                    auto rf = std::make_shared<RF>(cutoffLR, kappa, forceField, box, bc);
                    auto lj_rf = std::make_shared<LJ_RF>(kappa, forceField, box, bc, rf);
                    auto pair1 = std::make_pair(key1, lj_rf);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, lj_rf);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::RF) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    real_t kappa = properties::kappa(all);
                    auto rf = std::make_shared<RF>(cutoffLR, kappa, forceField, box, bc);
                    auto pair1 = std::make_pair(key1, rf);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, rf);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::LJ_SF) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto sf = std::make_shared<SF>(cutoffLR, forceField, box, bc);
                    auto lj_sf = std::make_shared<LJ_SF>(forceField, bc, sf);
                    auto pair1 = std::make_pair(key1, lj_sf);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, lj_sf);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::HS_SF) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto sf = std::make_shared<SF>(cutoffLR, forceField, box, bc);
                    auto hs_sf = std::make_shared<HS_SF>(forceField, bc, sf);
                    auto pair1 = std::make_pair(key1, hs_sf);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, hs_sf);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::HS_SC) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto screenedCoulomb = std::make_shared<SC>(forceField, bc);
                    auto hs_sc = std::make_shared<HS_SC>(forceField, bc, screenedCoulomb);
                    auto pair1 = std::make_pair(key1, hs_sc);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, hs_sc);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::HS_LEKNER) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto lekner = std::make_shared<Lekner>(box, bc);
                    auto hsLekner = std::make_shared<HardSphereLekner>(forceField, bc, lekner);
                    auto pair1 = std::make_pair(key1, hsLekner);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, hsLekner);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::SR) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto sr = std::make_shared<SoftRepulsion>(forceField, bc, cutoffSR);
                    auto pair1 = std::make_pair(key1, sr);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, sr);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::GAUSS_SF) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto gauss = std::make_shared<GaussianSF>(forceField, bc, cutoffLR, mesoscopic);
                    auto pair1 = std::make_pair(key1, gauss);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, gauss);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::GAUSS_SF_SR) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto gaussSF =
                            std::make_shared<GaussianSF>(forceField, bc, cutoffLR, mesoscopic);
                    auto sr =
                            std::make_shared<SoftRepulsion>(forceField, bc, cutoffSR);
                    auto gaussSF_SR =
                            std::make_shared<GaussianSF_SoftRepulsion>(forceField, bc, gaussSF, sr);
                    auto pair1 = std::make_pair(key1, gaussSF_SR);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, gaussSF_SR);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::NONE_INTERACTING) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto none = std::make_shared<NonInteracting>();
                    auto pair1 = std::make_pair(key1, none);
                    nonBondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, none);
                        nonBondedPairPotentials_.emplace(pair2);
                    }
                } else {
                    util::logAndThrow(logger,
                                      "Non-bonded pair potential '" + typeName +
                                      "' is not available.");
                }
            }

            // Bonded pair potentials.
            bondedPairPotentials_.clear();
            interactionParameters = forceField->bondedSpecifications();
            for (auto &ip: interactionParameters) {
                std::string typeName = ip.typeName;
                logger.debug(typeName + ": Bonded pair potential.");
                if (typeName == conf::HP) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto hp = std::make_shared<HP>(forceField, bc);
                    auto pair1 = std::make_pair(key1, hp);
                    bondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, hp);
                        bondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::HA_QP) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto ha_qp = std::make_shared<HalveAttractiveQP>(forceField, bc);
                    auto pair1 = std::make_pair(key1, ha_qp);
                    bondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, ha_qp);
                        bondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::HA_HP){
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                    auto ha_hp = std::make_shared<HalveAttractiveHP>(forceField, bc);
                    auto pair1 = std::make_pair(key1, ha_hp);
                    bondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, ha_hp);
                        bondedPairPotentials_.emplace(pair2);
                    }
                } else if (typeName == conf::HP_SR) {
                    std::string key1 = ip.spec1->name() + "-" + ip.spec2->name();
                     auto softRepulsion = std::make_shared<SoftRepulsion>(forceField, bc, cutoffSR);
                    auto harmonic = std::make_shared<HP>(forceField, bc);
                    auto hp_sr = std::make_shared<HarmonicSoftRepulsion>(forceField, bc, softRepulsion, harmonic);
                    auto pair1 = std::make_pair(key1, hp_sr);
                    bondedPairPotentials_.emplace(pair1);
                    if (ip.spec1 != ip.spec2) {
                        std::string key2 = ip.spec2->name() + "-" + ip.spec1->name();
                        auto pair2 = std::make_pair(key2, hp_sr);
                        bondedPairPotentials_.emplace(pair2);
                    }
                } else {
                    util::logAndThrow(logger, "Bonded pair potential '" + typeName + "' is not available.");
                }
            }

            // Log some info.
            auto nNonBonded = nonBondedPairPotentials_.size();
            logger.debug(std::to_string(nNonBonded) +
                         ": Number of associated non-bonded pair potentials.");
            auto nBonded = bondedPairPotentials_.size();
            logger.debug(std::to_string(nBonded) +
                         ": Number of associated bonded pair potentials.");

            logger.trace("Leaving.");
        }

        static void
        associateExternalPotentials(const ff_ptr_t &forceField,
                                    const box_ptr_t &box,
                                    const bc_ptr_t &bc,
                                    bool mesoscopic) {
            util::Logger logger("simploce::forces::associateExternalPotentials()");
            logger.trace("Entering.");

            auto external = forceField->externalSpecifications();
            for (const auto &spec: external) {
                std::string typeName = spec.typeName;
                if (typeName == conf::VOLTAGE) {
                    auto potential = std::make_shared<Voltage>(spec.e0, bc, spec.eps_r, mesoscopic);
                    externalPotentials_.emplace_back(potential);
                } else if (typeName == conf::WALL) {
                    auto plane = Plane::valueOf(spec.plane);
                    dist_t distance{spec.distance};
                    FlatSurface flatSurface{plane, distance};
                    srf_charge_density_t sigma{spec.sigma};
                    auto potential = std::make_shared<Wall>(spec.C12, spec.C6, bc, flatSurface, sigma);
                    externalPotentials_.emplace_back(potential);
                } else if (typeName == conf::PRESSURE_GRADIENT) {
                    auto potential = std::make_shared<PressureGradient>(spec.fe);
                    externalPotentials_.emplace_back(potential);
                } else if (typeName == conf::U_SRF_CG) {
                    Plane plane = Plane::valueOf(spec.plane);
                    FlatSurface flatSurface{plane};
                    auto potential = std::make_shared<UniformSurfaceChargeDensity>(spec.sigma,
                                                                                   flatSurface,
                                                                                   spec.eps_r,
                                                                                   bc,
                                                                                   spec.delta);
                    externalPotentials_.emplace_back(potential);
                } else if (typeName == conf::VIRTUAL_PLANES) {
                    auto potential = std::make_shared<VirtualPlanes>(box,
                                                                     bc,
                                                                     spec.spacing,
                                                                     spec.eps_r);
                    externalPotentials_.emplace_back(potential);
                } else {
                    util::logAndThrow(logger, "External potential of type'" + typeName + "' is not available.");
                }
            }

            // Log some info.
            logger.debug("Number of external potentials: " + std::to_string(externalPotentials_.size()));

            logger.trace("Leaving.");
        }

        /**
         * Returns non-bonded potential energy forces for given pairs of particles due to interactions
         * between these particles.
         * @return Potential energy and forces on particles.
         */
        static std::pair<energy_t, std::vector<force_t>>
        computeNonBondedForcesParticlePairs(const std::vector<PairList::p_pair_t> &particlePairs,
                                            const pp_map_t &nonBondedPairPotentials,
                                            int numberOfParticles) {
            static util::Logger logger("simploce::forces::computeNonBondedForcesParticlePairs()");
            logger.trace("Entering.");

            static std::size_t counter = 0;

            // Compute interactions for all particle pairs.
            std::vector<force_t> forces(numberOfParticles, force_t{0.0, 0.0, 0.0});
            energy_t energy{0.0};

            if (particlePairs.empty()) {
                return std::move(std::make_pair(energy, forces));
            }

            for (auto &particlePair: particlePairs) {
                // First particle
                p_ptr_t pi = particlePair.first;
                auto index_i = pi->index();

                // Second particle
                p_ptr_t pj = particlePair.second;
                auto index_j = pj->index();

                // Call pair potential.
                std::string key = pi->spec()->name() + "-" + pj->spec()->name();
                pair_potential &pairPotential = *selectPairPotential(key, nonBondedPairPotentials);
                std::pair<energy_t, force_t> ef = pairPotential(pi, pj);

                // Store results.
                energy += std::get<0>(ef);
                forces[index_i] += std::get<1>(ef);
                forces[index_j] -= std::get<1>(ef);

                // Next.
            }

            // Done.
            counter += 1;

            logger.trace("Leaving.");
            return std::make_pair(energy, forces);
        }

        /**
         * Computes non-bonded forces and associated interaction energies.
         * @return Non-bonded potential energy.
         */
        static energy_t
        computeNonBondedForces(bool concurrent,
                               const std::vector<p_ptr_t> &all,
                               const pairlist_ptr_t &pairList,
                               const pp_map_t &nonBondedPairPotentials) {

            static util::Logger logger("simploce::forces::computeNonBondedForces()");
            logger.trace("Entering.");

            static std::vector<std::vector<PairList::p_pair_t>> subParticlePairs;
            static bool firstTime = true;

            if (pairList->isEmpty()) {
                return 0.0;
            }

            auto numberOfParticles = pairList->numberOfParticles();
            std::vector<result_t> results{};

            if (numberOfParticles > simploce::conf::MIN_NUMBER_OF_PARTICLES && concurrent) {
                // Handle particle-particle interaction concurrently as tasks,
                // where one task is executed by the current thread.
                if (firstTime) {
                    logger.info("Non-bonded forces and energies are computed concurrently.");
                }
                std::vector<std::future<result_t> > futures{};
                if (pairList->isModified() || firstTime) {
                    subParticlePairs = util::makeSubLists(pairList->nonBondedParticlePairs());
                }
                auto numberOfTasks = subParticlePairs.size() - 1;
                if (numberOfTasks > 1) {
                    // Set up concurrent force calculations.
                    for (std::size_t k = 0; k != numberOfTasks; ++k) {
                        const auto &single = subParticlePairs[k];
                        futures.push_back(
                                std::async(
                                        std::launch::async,
                                        computeNonBondedForcesParticlePairs,
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
                const auto &single = *(subParticlePairs.end() - 1);
                if (!single.empty()) {
                    auto result = computeNonBondedForcesParticlePairs(single,
                                                                      nonBondedPairPotentials,
                                                                      int(numberOfParticles));
                    results.push_back(result);
                }
            } else {
                // Sequentially. All particle-particle pair interactions are handled a single thread.
                if (firstTime) {
                    logger.info("Non-bonded forces and energies are computed sequentially.");
                }
                auto result =
                        computeNonBondedForcesParticlePairs(pairList->nonBondedParticlePairs(),
                                                            nonBondedPairPotentials,
                                                            int(numberOfParticles));
                results.push_back(result);
            }

            // Collect potential energies and forces. Update forces on particles.
            energy_t energy{0.0};
            for (auto &result: results) {
                const auto &forces = result.second;
                for (auto &particle: all) {
                    auto index = particle->index();
                    force_t f = forces[index] + particle->force();
                    particle->force(f);
                }
                energy += result.first;
            }

            // Done. Returns potential energy.
            firstTime = false;

            logger.trace("Leaving.");
            return energy;
        }


        /**
         * Computes bonded forces and associated potential energies. Forces are added to particles.
         * @return Bonded interaction potential energy.
         */
        static energy_t
        computeBondedForces(const std::vector<pg_ptr_t> &groups,
                            const pp_map_t &bondedPairPotentials) {
            static util::Logger logger{"simploce::forces::computeBondedForces()"};
            logger.trace("Entering");

            energy_t bonded{0.0};
            for (auto &group: groups) {
                for (auto &bond: group->bonds()) {
                    // Get bonded particles.
                    auto p1 = bond.getParticleOne();
                    auto p2 = bond.getParticleTwo();

                    // Call pair potential.
                    std::string key = p1->spec()->name() + "-" + p2->spec()->name();
                    pair_potential &pairPotential = *selectPairPotential(key, bondedPairPotentials);
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

            logger.trace("Leaving.");
            return bonded;
        }

        static energy_t
        computeExternalForces(const std::vector<p_ptr_t> &all,
                              const ep_cont_t &externalPotentials) {
            static util::Logger logger{"simploce::forces::computeExternalForces()"};
            logger.trace("Entering.");

            energy_t external{};
            for (auto &particle: all) {
                for (auto &ep: externalPotentials) {
                    auto result = (*ep)(particle);
                    external += result.first;
                    auto f = particle->force();
                    f += result.second;
                    particle->force(f);
                }
            }

            logger.trace("Leaving.");
            return external;
        }

        static void
        initiateExternalPotentials(const p_system_ptr_t& particleSystem,
                                   const ep_cont_t &externalPotentials) {
            static util::Logger logger("simploce::forces::initiateExternalPotentials");
            logger.trace("Entering");

            for (auto& ep: externalPotentials) {
                ep->initialize(particleSystem);
            }

            logger.trace("Leaving");
        }

        /**
         * Update external potentials, whatever is required for a proper
         * functioning of these potentials.
         * @param particleSystem Particle system
         * @param externalPotentials External potentials.
         */
        static void
        updateExternalPotentials(const p_system_ptr_t& particleSystem,
                                 const ep_cont_t &externalPotentials) {
            static util::Logger logger("simploce::forces::updateExternalPotentials");
            logger.trace("Entering");

            for (auto& ep: externalPotentials) {
                ep->update(particleSystem);
            }

            logger.trace("Leaving");
        }

        /**
         * Update external potentials, whatever is required for a proper
         * functioning of these potentials.
         * @param particle Particle
         * @param externalPotentials External potentials.
         */
        static void
        updateExternalPotentials(const p_ptr_t& particle,
                                 const ep_cont_t &externalPotentials) {
            for (auto& ep: externalPotentials) {
                ep->update(particle);
            }
        }

        /**
         * Computer interaction energy of given particle with other particles.
         * @param particle Given particle.
         * @param other Other particles. Given particle should not be in this collection.
         * @param nonBondedPairPotentials Non-bonded potential.
         * @param cutoffs Cutoff distances.
         * @param bc Boundary condition.
         * @return Non-bonded interaction energy.
         */
        static energy_t
        computeNonBoundedInteractionEnergyOfOneParticle(const p_ptr_t &particle,
                                                        const std::vector<p_ptr_t> &other,
                                                        const pp_map_t &nonBondedPairPotentials,
                                                        const rc_ptr_t &cutoffs,
                                                        const bc_ptr_t &bc) {
            static util::Logger logger{"simploce::forces::computeNonBoundedInteractionEnergyOfOneParticle()"};
            logger.trace("Entering.");

            static auto cutoff = cutoffs->longRanged();
            static auto rc2 = cutoff() * cutoff();

            energy_t nonBonded{0.0};
            auto name = particle->spec()->name();
            auto id = particle->id();
            auto pr = particle->position();

            // Non-bonded interaction with free particles.
            for (const auto &p: other) {
                if (p->id() != id) {  // No self interaction.
                    auto r = p->position();
                    dist_vect_t R = bc->apply(pr, r);
                    auto R2 = norm_square<real_t>(R);
                    if (R2 <= rc2) {
                        std::string key = name + "-" + p->spec()->name();
                        auto &pairPotential = *nonBondedPairPotentials.at(key);

                        // Force is ignored.
                        std::pair<energy_t, force_t> ef = pairPotential(particle, p);
                        nonBonded += std::get<0>(ef);
                    }
                }
            }
            logger.trace("Leaving");
            return nonBonded;
        }

        /**
         * Returns interaction energy of a given particle with all other particles
         * (except with itself) within the cutoff distance. The cutoff distance corresponds
         * to the cutoff distance for long-ranged interactions.
         * @param particle Particle.
         * @param all -All- particles.
         * @param nonBondedPairPotentials Pair potentials.
         * @return Bonded, non-bonded, and external potential energy.
         */
        static std::tuple<energy_t, energy_t, energy_t>
        interactionEnergyOfOneParticle(const p_ptr_t &particle,
                                       const std::vector<p_ptr_t> &free,
                                       const std::vector<pg_ptr_t> &groups,
                                       const pp_map_t &nonBondedPairPotentials,
                                       const pp_map_t &bondedPairPotentials,
                                       const ep_cont_t &externalPotentials,
                                       const rc_ptr_t &cutoffs,
                                       const bc_ptr_t &bc) {
            static util::Logger logger{"simploce::forces::interactionEnergyOfOneParticle()"};
            logger.trace("Entering.");

            static auto cutoff = cutoffs->longRanged();
            static auto rc2 = cutoff() * cutoff();
            static bool firstTime{true};
            static std::vector<std::vector<p_ptr_t>> particleSets;
            static size_t numberOfTasks = 0;

            energy_t nonBonded{0.0}, bonded{0.0};
            auto name = particle->spec()->name();
            auto id = particle->id();
            auto pr = particle->position();

            // Non-bonded interaction with free particles.
            logger.trace("Calculating non-bonded interaction energy.");
            std::vector<energy_t> results{};
            if (free.size() > simploce::conf::MIN_NUMBER_OF_PARTICLES) {
                if (firstTime) {
                    logger.info("Non-bonded interaction energy is computed concurrently.");
                    particleSets = util::makeSubLists(free);
                    numberOfTasks = particleSets.size() - 1;
                }
                std::vector<std::future<energy_t>> futures{};
                if (numberOfTasks > 1) {
                    for (auto k = 0; k != numberOfTasks; ++k) {
                        const auto& particleSet = particleSets[k];
                        futures.push_back(std::async(
                            std::launch::async,
                            computeNonBoundedInteractionEnergyOfOneParticle,
                            std::ref(particle),
                            std::ref(particleSet),
                            std::ref(nonBondedPairPotentials),
                            std::ref(cutoffs),
                            std::ref(bc)
                        ));
                    }
                    // Wait for each task to complete.
                    results = util::waitForAll(futures);
                }
                // One remaining set of particleSet is handled by the current thread.
                const auto& particleSet = particleSets[particleSets.size() - 1];
                if (!particleSet.empty()) {
                    auto energy = computeNonBoundedInteractionEnergyOfOneParticle(particle,
                                                                                  particleSet,
                                                                                  nonBondedPairPotentials,
                                                                                  cutoffs,
                                                                                  bc);
                    results.emplace_back(energy);
                }
            } else {
                if (firstTime) {
                    logger.info("Non-bonded interaction energy is computed sequentially.");
                }
                auto energy = computeNonBoundedInteractionEnergyOfOneParticle(particle,
                                                                              free,
                                                                              nonBondedPairPotentials,
                                                                              cutoffs,
                                                                              bc);
                results.emplace_back(energy);
            }
            // Collect the results.
            for (auto& r: results) {
                nonBonded += r;
            }

            /*
            for (const auto &p: free) {
                if (p->id() != id) {
                    auto r = p->position();
                    dist_vect_t R = bc->apply(pr, r);
                    auto R2 = norm_square<real_t>(R);
                    if (R2 <= rc2) {
                        std::string key = name + "-" + p->spec()->name();
                        auto &pairPotential = *nonBondedPairPotentials.at(key);
                        // Force is ignored.
                        std::pair<energy_t, force_t> ef = pairPotential(particle, p);
                        nonBonded += std::get<0>(ef);
                    }
                }
            }
            */

            // Interaction with other and within groups.
            for (const auto &g: groups) {
                if (!g->contains(particle)) {
                    auto particles = g->particles();
                    for (const auto &p: particles) {
                        auto r = p->position();
                        dist_vect_t R = bc->apply(pr, r);
                        auto R2 = norm_square<real_t>(R);
                        if (R2 <= rc2) {
                            std::string key = name + "-" + p->spec()->name();
                            auto &pairPotential = *nonBondedPairPotentials.at(key);
                            // Force is ignored.
                            std::pair<energy_t, force_t> ef = pairPotential(particle, p);
                            nonBonded += std::get<0>(ef);
                        }
                    }
                } else {
                    auto particles = g->particles();
                    for (const auto &p: particles) {
                        if (particle->id() != p->id()) {
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
            logger.trace("Calculating energy due to external potentials.");
            energy_t external{0.0};
            for (auto &ep: externalPotentials) {
                auto &potential = *ep;
                std::pair<energy_t, force_t> ef = potential(particle);
                // Force is ignored.
                external += std::get<0>(ef);
            }

            logger.trace("Leaving.");
            firstTime = false;
            return std::move(std::make_tuple(bonded, nonBonded, external));
        }

        static void
        associatePotentials(const rc_ptr_t &cutoffs,
                            const std::vector<p_ptr_t> &all,
                            const ff_ptr_t &forceField,
                            bool mesoscopic,
                            const box_ptr_t &box,
                            const bc_ptr_t &bc) {
            static bool ASSOCIATED{false};
            if (!ASSOCIATED) {
                associatePairPotentials(cutoffs, all, forceField, mesoscopic, box, bc);
                associateExternalPotentials(forceField, box, bc, mesoscopic);
                ASSOCIATED = true;
            }
        }

        static void
        associatePotentials(const p_system_ptr_t& particleSystem,
                            const rc_ptr_t& cutoffs,
                            const bc_ptr_t& bc,
                            const ff_ptr_t& forceField,
                            bool mesoscopic) {
            static bool ASSOCIATED{false};
            if ( !ASSOCIATED ) {
                auto box = particleSystem->box();
                particleSystem->doWithAll<void>([cutoffs,
                                                      box,
                                                      bc,
                                                      forceField,
                                                      mesoscopic] (const std::vector<p_ptr_t>& all) {
                    forces::associatePotentials(cutoffs,
                                                all,
                                                forceField,
                                                mesoscopic,
                                                box,
                                                bc);
                    ASSOCIATED = true;
                });
            }
        }
    }

    Forces::Forces(param_ptr_t param,
                   bc_ptr_t bc,
                   ff_ptr_t forceField) :
            param_{std::move(param)}, cutoffs_{}, bc_{std::move(bc)},
            forceField_{std::move(forceField)} {
        util::Logger logger{"simploce::Forces::Forces()"};
        logger.trace("Entering.");

        if (!param_) {
            logAndThrow(logger, "Missing parameters.");
        }
        auto ignore = param_->get<bool>("simulation.forces.ignore_cutoff", false);
        if (!ignore) {
            dist_t sr = param_->get<real_t>("simulation.forces.cutoffSR");
            dist_t lr = param_->get<real_t>("simulation.forces.cutoffLR");
            cutoffs_ = std::make_shared<Cutoffs>(sr, lr);
        } else {
            cutoffs_ = std::make_shared<Cutoffs>(conf::LARGE, conf::LARGE);
            logger.warn("Cutoff distances are ignored (all set to a large value).");
        }
        concurrent_ = param_->get<bool>("simulation.forces.concurrent", true);
        mesoscopic_ = param_->get<bool>("simulation.mesoscale");
        if (!cutoffs_) {
            logAndThrow(logger, "Missing cutoff distances.");
        }
        if (!bc_) {
            logAndThrow(logger, "Missing boundary conditions.");
        }
        if (!forceField_) {
            logAndThrow(logger, "Missing force field.");
        }

        logger.trace("Leaving.");
    }

    energy_t
    Forces::nonBonded(const p_system_ptr_t& particleSystem,
                      const pairlist_ptr_t& pairList) {
        static util::Logger logger("simploce::Forces::nonBonded()");
        logger.trace("Entering.");

        // Must be done first, if not done already.
        forces::associatePotentials(particleSystem,
                                    this->cutoffs_,
                                    bc_,
                                    forceField_,
                                    mesoscopic_);

        static bool firstTime = true;
        if (firstTime) {
            logger.debug(std::to_string(cutoffs_->shortRanged()()) +
            ": Cutoff distance for short ranged interactions.");
            logger.debug(std::to_string(cutoffs_->longRanged()()) +
            ": Cutoff distance for long ranged interactions.");
            firstTime = false;
        }

        auto box = particleSystem->box();
        auto energy =
            particleSystem->doWithAll<energy_t>([this, pairList, box] (std::vector<p_ptr_t>& all) {
                forces::associatePotentials(cutoffs_, all, forceField_, mesoscopic_, box, this->bc_);
                return forces::computeNonBondedForces(this->concurrent_,
                                                      all,
                                                      pairList,
                                                      forces::nonBondedPairPotentials_);
            });

        logger.trace("Leaving");
        return energy;
    }

    energy_t
    Forces::external(const p_system_ptr_t& particleSystem) {
        static util::Logger logger{"simploce::Forces::external()"};
        logger.trace("Entering.");

        // Must be done first, if not done already.
        forces::associatePotentials(particleSystem,
                                    this->cutoffs_,
                                    bc_,
                                    forceField_,
                                    mesoscopic_);

        auto box = particleSystem->box();
        auto energy = particleSystem->doWithAll<energy_t>([this, box] (const std::vector<p_ptr_t>& all) {
            forces::associatePotentials(this->cutoffs_,
                                        all, this->forceField_,
                                        this->mesoscopic_, box,
                                        this->bc_);
            return forces::computeExternalForces(all, forces::externalPotentials_);
        });

        logger.trace("Leaving");
        return energy;
    }

    energy_t
    Forces::bonded(const p_system_ptr_t& particleSystem) {
        // Must be done first, if not done already.
        forces::associatePotentials(particleSystem,
                                    this->cutoffs_,
                                    bc_,
                                    forceField_,
                                    mesoscopic_);

        auto box = particleSystem->box();
        return particleSystem->doWithAllFreeGroups<energy_t>([this, box] (
                const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups
                ) {
            forces::associatePotentials(this->cutoffs_,
                                        all,
                                        this->forceField_,
                                        this->mesoscopic_,
                                        box,
                                        this->bc_);
            return forces::computeBondedForces(groups,
                                               forces::bondedPairPotentials_);
        });
    }

    std::tuple<energy_t, energy_t, energy_t>
    Forces::interaction(const p_ptr_t& particle,
                        const p_system_ptr_t& particleSystem) {
        static util::Logger logger("simploce::Forces::interaction()");
        logger.trace("Entering");

        // Must be done first, if not done already.
        forces::associatePotentials(particleSystem,
                                    this->cutoffs_,
                                    bc_,
                                    forceField_,
                                    mesoscopic_);

        // Calculate interaction energy.
        auto box = particleSystem->box();
        auto result =
                particleSystem->doWithAllFreeGroups<std::tuple<energy_t, energy_t, energy_t>>(
                [this, particle, box] (
                    const std::vector<p_ptr_t>& all,
                    const std::vector<p_ptr_t>& free,
                    const std::vector<pg_ptr_t>& groups) {
            return forces::interactionEnergyOfOneParticle(particle,
                                                          free,
                                                          groups,
                                                          forces::nonBondedPairPotentials_,
                                                          forces::bondedPairPotentials_,
                                                          forces::externalPotentials_,
                                                          this->cutoffs_,
                                                          this->bc_);
        });

        logger.trace("Leaving");
        return std::move(result);
    }

    void
    Forces::initiate(const simploce::p_system_ptr_t &particleSystem) {
        // Must be done first, if not done already.
        forces::associatePotentials(particleSystem,
                                    cutoffs_,
                                    bc_,
                                    forceField_,
                                    mesoscopic_);
        forces::initiateExternalPotentials(particleSystem,
                                           forces::externalPotentials_);
    }

    void
    Forces::update(const simploce::p_system_ptr_t &particleSystem) {
        forces::updateExternalPotentials(particleSystem,
                                         forces::externalPotentials_);
    }

    void
    Forces::update(const simploce::p_ptr_t &particle) {
         forces::updateExternalPotentials(particle,
                                         forces::externalPotentials_);
    }

    void
    Forces::fallback() {
        for (auto& ep: forces::externalPotentials_) {
            ep->fallback();
        }
    }

    void
    Forces::complete() {
        for (auto& ep: forces::externalPotentials_) {
            ep->complete();
        }
    }

}
