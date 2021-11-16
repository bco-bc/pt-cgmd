/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on on 11/11/21.
 */

#include "simploce/simulation/lj-electrostatics.hpp"
#include "simploce/simulation/force-field.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/util/util.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/types/cvector_t.hpp"
#include <vector>
#include <future>
#include <memory>
#include <tuple>
#include <cassert>
#include <utility>

namespace simploce {
    namespace forces {

        using result_t = std::pair<energy_t, std::vector<force_t>>;

        /**
         * Returns interaction potential energy and force on particle i. The Coulomb
         * interaction between particle i and j is calculated according to the
         * Shifted Force (SF) method of
         * Levitt, M. et al, Comput. Phys. Commun. 1995, 91, 215−231.
         * @param ri Position particle i.
         * @param qi Charge particle j.
         * @param rj Position particle j.
         * @param qj Charge particle j.
         * @param C12 Lennard-Jones C12 parameter.
         * @param C6 Lennard-Jones C6 parameters.
         * @param eps_r Relative permittivity.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Potential energy, force on particle i, and the applied distance between i and j.
         */
        static std::tuple<energy_t, force_t, length_t>
        forceSF_(const position_t& ri,
                 const charge_t& qi,
                 const position_t& rj,
                 const charge_t& qj,
                 real_t C12,
                 real_t C6,
                 real_t eps_r,
                 const box_ptr_t &box,
                 const bc_ptr_t& bc) {
            // Initializations.
            static const real_t four_pi_e0 = units::mu<real_t>::FOUR_PI_E0;
            static const length_t rc = properties::cutoffDistance(box);
            static const real_t rc2 = rc() * rc();

            // Apply boundary condition.
            dist_vect_t rij = bc->apply(ri, rj);
            auto Rij = norm<real_t>(rij);

            // Potential energy
            // LJ
            real_t Rij2 = Rij * Rij;
            real_t Rij6 = Rij2 * Rij2 * Rij2;
            real_t Rij12 = Rij6 * Rij6;
            real_t t1 = C12 / Rij12;
            real_t t2 = C6 / Rij6;
            real_t LJ = t1 - t2;  // kj/mol
            // SF
            real_t t3 = qi * qj / (four_pi_e0 * eps_r);
            real_t SF = t3 * (1.0 / Rij - 1.0 / rc() + (Rij - rc()) / rc2);  // kJ/mol
            // Total potential energy.
            energy_t epot = LJ + SF;

            // Forces.
            dist_vect_t uv = rij/Rij;                          // Unit vector
            real_t dLJdR = -6.0 * ( 2.0 * t1 - t2 ) / Rij;     // LJ, kJ/(mol nm).
            real_t dSFdR = t3 * (-1.0 / Rij2 + 1.0 / rc2);     // SF, kJ/(mol nm)
            force_t f{};
            for (std::size_t k = 0; k != 3; ++k) {
                real_t fLJ = -dLJdR * uv[k];                   // kJ/(mol nm)
                real_t fElec = -dSFdR * uv[k];                 // kJ/(mol nm)
                f[k] = fLJ + fElec;
            }
            // Done
            return std::make_tuple(epot, f, Rij);
        }

        template <typename P>
        static std::pair<energy_t, std::vector<force_t>>
        ppForces_(const typename PairLists<P>::pp_pair_cont_t &particlePairs,
                  const ff_ptr_t &forceField,
                  int numberOfParticles,
                  const box_ptr_t &box,
                  const bc_ptr_t &bc) {
            // Particle pointer type.
            using p_ptr_t = std::shared_ptr<P>;

            // Interaction parameters.
            auto ljParameters = forceField->lennardJonesParameters();
            auto eps_r = forceField->relativePermittivity();

            // Compute interactions for all particle pairs.
            std::vector<force_t> forces(numberOfParticles, force_t{});
            energy_t epot{0.0};
            for (auto &particlePair : particlePairs) {

                // First particle
                p_ptr_t pi = particlePair.first;
                position_t ri = pi->position();
                charge_t qi = pi->charge();
                auto index_i = pi->index();

                // Second particle
                p_ptr_t pj = particlePair.second;
                position_t rj = pj->position();
                charge_t qj = pj->charge();
                auto index_j = pj->index();

                // Interaction parameters
                auto param = ljParameters.at(pi->spec(), pj->spec());
                auto C12 = param.first;
                auto C6 = param.second;

                // Replaceable.
                // Lennard-Jones + Shifted Force electrostatic interaction.
                auto ef = forceSF_(ri, qi, rj, qj, C12, C6, eps_r, box, bc);

#ifdef _DEBUG
                // Too close?
                util::tooClose<P>(pi, pj, ef);
#endif
                // Store energy and forces.
                epot += std::get<0>(ef);
                forces[index_i] += std::get<1>(ef);
                forces[index_j] -= std::get<1>(ef);
            }
            // Done.
            return std::make_pair(epot, forces);
        }

        template <typename P>
        static energy_t
        compute_(const std::vector<std::shared_ptr<P>> &all,
                 const PairLists<P> &pairLists,
                 const ff_ptr_t &forceField,
                 const box_ptr_t &box,
                 const bc_ptr_t &bc) {

            using pp_pair_cont_t = typename PairLists<P>::pp_pair_cont_t;

            static util::Logger logger("simploce::forces::compute_()");
            static std::vector<pp_pair_cont_t> subPairLists;
            static bool firstTime = true;

            auto numberOfParticles = all.size();
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
                                ppForces_<P>,
                                std::ref(single),
                                std::ref(forceField),
                                numberOfParticles,
                                std::ref(box),
                                std::ref(bc)
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
                    auto result = ppForces_<P>(single, forceField, numberOfParticles, box, bc);
                    results.push_back(result);
                }
            } else {
                // Sequentially. All particle-particle pair interactions are handled by the
                // current thread.
                if ( firstTime ) {
                    logger.debug("Lennard-Jones and electrostatic forces and energies are computed concurrently.");
                }
                auto result = ppForces_<P>(pairLists.particlePairList(),
                                           forceField,
                                           numberOfParticles,
                                           box,
                                           bc);
                results.push_back(result);
            }

            // Collect potential energies and forces.
            energy_t epot{0.0};
            for (auto& result : results) {
                const auto& forces = result.second;
                for (auto& particle : all) {
                    auto index = particle->index();
                    force_t f = forces[index] + particle->force();
                    particle->force(f);
                }
                epot += result.first;
            }

            // Done. Returns potential energy.
            firstTime = false;
            return epot;
        }

        energy_t computeLjElectrostatic(const std::vector<atom_ptr_t>& all,
                                        const PairLists<Atom> &pairLists,
                                        const ff_ptr_t& forceField,
                                        const box_ptr_t &box,
                                        const bc_ptr_t &bc) {
            return std::move(compute_<Atom>(all, pairLists, forceField, box, bc));
        }

        energy_t computeLjElectrostatic(const std::vector<bead_ptr_t>& all,
                                        const PairLists<Bead> &pairLists,
                                        const ff_ptr_t& forceField,
                                        const box_ptr_t &box,
                                        const bc_ptr_t &bc) {
            return std::move(compute_<Bead>(all, pairLists, forceField, box, bc));
        }
    }
}

