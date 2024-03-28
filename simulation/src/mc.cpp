/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 19 October 2019, 18:56
 */
#include "simploce/simulation/mc.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/util.hpp"
#include "simploce/units/units-mu.hpp"
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <random>
#include <utility>
#include <cmath>
#include <tuple>
#include <vector>

namespace simploce {
    namespace mc {

        static const real_t LIMIT = 75.0;

        /**
         * Returns coordinate descriptor.
         * @param k Index, must be 0 <= k< =2.
         * @return "x", or "y" or "z".
         */
        static std::string
        coordinateDescriptor_(std::size_t k) {
            if (k > 2) {
                throw std::domain_error(std::to_string(k) + ": Illegal value for coordinate descriptor");
            }
            std::string c;
            if (k == 0) c = "x";
            if (k == 1) c = "y";
            if (k == 2) c = "z";
            return c;
        }

        /**
         * Returns coordinate index.
         * @param c Coordinate descriptor.
         * @return Index (0, 1, or 2).
         */
        static int
        coordinateIndex_(const std::string& c) {
            if (c == "x") return 0;
            if (c == "y") return 1;
            if (c == "z") return 2;
            throw std::domain_error(c + ": Illegal value for coordinate index.");
        }

        /**
         * Returns total potential energy.
         * @param result Current values for non-bonded, bonded, and external interactions.
         * @return Sum of potential energies.
         */
        static energy_t
        sumPotentialEnergies_(const std::tuple<energy_t, energy_t, energy_t> &result) {
            return std::get<0>(result) + std::get<1>(result) + std::get<2>(result);
        }

        /**
         * Accept the current state.
         * @param result_i Initial state.
         * @param result_f Final state.
         * @return Result.
         */
        static std::tuple<energy_t, energy_t, energy_t, bool>
        accept_(const std::tuple<energy_t, energy_t, energy_t> &result_i,
                const std::tuple<energy_t, energy_t, energy_t> &result_f) {
            return std::make_tuple(
                    std::get<0>(result_f) - std::get<0>(result_i),
                    std::get<1>(result_f) - std::get<1>(result_i),
                    std::get<2>(result_f) - std::get<2>(result_i),
                    true
            );
        }

        /**
         * Finds three booleans associated with three coordinates whether to place or keep a
         * given x, y or z-coordinate or combinations of these in the central box.
         * @param s string of the form "x,y,z", "x,y", "y,z", "z", etc.
         * @return Vector of 3 booleans.
         */
        static std::vector<bool>
        keepInsideBox_(std::string s) {
            util::Logger logger("simulation::mc::keepInsideBox_");
            logger.trace("Entering.");
            logger.debug(s + ": value of coordinate string");
            std::vector<bool> inside(3, false);
            if (s.empty()) {
                return inside;
            }
            boost::trim(s);
            std::vector<std::string> results;
            boost::split(results, s, boost::is_any_of(","));
            for (auto & r : results) {
                auto j = coordinateIndex_(r);
                inside[j] = true;
            }
            logger.trace("Leaving.");
            return inside;
        }

        /**
         * Displaces the given particle.
         * @return Difference in bonded, non-bonded potential and external potential energy, and
         * acceptance.
         */
        static std::tuple<energy_t, energy_t, energy_t, bool>
        displaceParticle_(p_ptr_t &particle,
                          const p_system_ptr_t &particleSystem,
                          const interactor_ptr_t &interactor,
                          const param_ptr_t &param) {
            static util::Logger logger{"simploce::mc::displaceParticle_"};
            logger.trace("Entering");

            // Setup
            static int counter = 0;
            static temperature_t temperature = param->get<real_t>("simulation.temperature");
            static auto range = param->get<real_t>("simulation.displacer.mc.range");
            static auto mesoscopic = param->get<bool>("simulation.mesoscale");
            static auto inBox = param->get<bool>("simulation.displacer.mc.in-box", false);
            static auto inside = keepInsideBox_(
                param->get<std::string>("simulation.displacer.mc.keep-in-box", "")
            );
            static auto zNonNegative = param->get<bool>("simulation.displacer.mc.z-non-negative", false);
            static auto KB = mesoscopic ? 1.0 : units::mu<real_t>::KB;
            static const real_t kT = KB * temperature();
            static std::vector<real_t> ranges(3);

            if (counter == 0) {
                logger.info(
                    std::to_string(inBox) +
                    ": Always put selected particle in the simulation box?"
                );
                for (std::size_t k = 0; k != 3; ++k) {
                    auto c = coordinateDescriptor_(k);
                    logger.info(std::to_string(inside[k]) +
                                ": Keep selected particle's " + c + "-coordinate in simulation box?");
                }
                logger.info(std::to_string(zNonNegative) + ": z-coordinate must be a non-negative number?");

                // Find ranges for displacement.
                auto box = *particleSystem->box();
                logger.debug(util::to_string(box) + ": Box dimensions");
                for (std::size_t k = 0; k != 3; ++k) {
                    ranges[k] = inBox || inside[k] ? box[k] : range;
                    auto c = coordinateDescriptor_(k);
                    logger.debug(std::to_string(ranges[k]) + ": Displacement range for " + c + "-coordinate.");
                }
            }
            counter += 1;

            // Random number generators.
            std::random_device rd;
            std::mt19937_64 gen(rd());
            std::uniform_real_distribution<real_t> dis01(0.0, 1.0);

            // Current position.
            auto ri = particle->position();
            auto pri = particle->previousPosition();  // Previous position.

            // Calculate initial (current) energy.
            auto result_i = interactor->interact(particle, particleSystem);
            auto energy_i = sumPotentialEnergies_(result_i);
            if (energy_i() >= conf::LARGE || std::isnan(energy_i())) {
                energy_i = conf::LARGE;
            }

            // Displace/move particle. Update state.
            position_t rf; // Displaced/final position.
            for (std::size_t k = 0; k != 3; ++k) {
                std::uniform_real_distribution<real_t> disCoordinate(0.0, ranges[k]);
                auto rf_k = inBox || inside[k] ?
                        disCoordinate(gen) :
                        ri[k] - 0.5 * ranges[k] + disCoordinate(gen);
                if (k == 2 && zNonNegative) {
                    rf_k = rf_k < 0.0 ? -rf_k : rf_k;
                }
                rf[k] = rf_k;
            }
            particle->position(rf);             // Sets new position.
            particle->previousPosition(ri);     // Sets previous position.
            interactor->update(particle);

            // Calculate final (new) energy.
            auto result_f = interactor->interact(particle, particleSystem);
            auto energy_f = sumPotentialEnergies_(result_f);
            if (energy_f() >= conf::LARGE || std::isnan(energy_f())) {
                energy_f = 2.0 * conf::LARGE;
            }

            // Accept or reject displacement.
            energy_t difference = energy_f() - energy_i();
            real_t difference_over_kT = difference() / kT;
            if (difference_over_kT < LIMIT && energy_f() < conf::LARGE) {
                if (difference_over_kT > 0.0) {
                    real_t w = std::exp(-difference_over_kT);
                    real_t rv = dis01(gen);
                    if (rv > w) {
                        // Reject. Restore previous state.
                        particle->position(ri);
                        particle->previousPosition(pri);
                        interactor->fallback();
                        return std::make_tuple(0.0, 0.0, 0.0, false);
                    } else {
                        // Accept. Keep new state.
                        return std::move(accept_(result_i, result_f));
                    }
                } else {
                    // Accept. Keep new state.
                    return std::move(accept_(result_i, result_f));
                }
            } else {
                // Reject. Restore original state.
                particle->position(ri);
                particle->previousPosition(pri);
                interactor->fallback();
                return std::make_tuple(0.0, 0.0, 0.0, false);
            }
        }

        static SimulationData
        displaceOneParticle_(const std::vector<p_ptr_t> &displaceables,
                             const p_system_ptr_t &particleSystem,
                             const interactor_ptr_t &interactor,
                             const param_ptr_t &param) {
            static util::Logger logger("simploce::displaceOneParticle_");
            logger.trace("Entering.");
            static int counter = 0;

            if (counter == 0) {
                logger.info(std::to_string(displaceables.size()) +
                            ": Number of displaceable particles.");
                counter += 1;
            }

            // Set up.
            static temperature_t temperature = param->get<real_t>("simulation.temperature");
            static energy_t bonded = 0.0;
            static energy_t nonBonded = 0.0;
            static energy_t external = 0.0;

            // Randomly select one particle.
            std::random_device rd;
            std::mt19937_64 gen(rd());
            std::uniform_int_distribution<std::size_t> dis(0, displaceables.size() - 1);
            auto index = dis(gen);
            auto particle = displaceables[index];
            logger.debug(util::to_string(particle->id()) +
                         ": Identifier selected particle.");

            // Displace selected particle.
            auto result =
                    displaceParticle_(particle, particleSystem, interactor, param);

            // Store some data.
            SimulationData data;
            bonded += std::get<0>(result);
            nonBonded += std::get<1>(result);
            external += std::get<2>(result);
            data.bonded = bonded;
            data.nonBonded = nonBonded;
            data.external = external;
            data.accepted = std::get<3>(result);
            data.temperature = temperature;

            logger.trace("Leaving.");
            return data;
        }

    }
    
    MonteCarlo::MonteCarlo(param_ptr_t param,
                           interactor_ptr_t interactor) :
            param_{std::move(param)}, interactor_{std::move(interactor)} {
    }

    SimulationData
    MonteCarlo::displace(const p_system_ptr_t& particleSystem) const {
        static std::vector<p_ptr_t> displaceables;
        return particleSystem->doWithAll<SimulationData>([this, particleSystem] (
                const std::vector<p_ptr_t>& particles) {
            if (particleSystem->changed() || displaceables.empty()) {
                displaceables.clear();
                for (auto &p: particles) {
                    if (!p->frozen())
                        displaceables.emplace_back(p);
                }
            }
            return mc::displaceOneParticle_(displaceables,
                                            particleSystem,
                                            this->interactor_,
                                            this->param_);
        });
    }

}