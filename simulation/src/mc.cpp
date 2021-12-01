/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 19 October 2019, 18:56
 */
#include "simploce/simulation/mc.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/util.hpp"
#include "simploce/units/units-mu.hpp"
#include <random>
#include <utility>
#include <cmath>
#include <tuple>

namespace simploce {
    
    static const real_t RANGE = 0.1;
    static const real_t LIMIT = 75.0;
    
    // Returns difference in bonded and non-bonded potential energy, and acceptance.
    static std::tuple<energy_t, energy_t, bool> 
    displaceParticle_(p_ptr_t& particle,
                      const p_system_ptr_t& particleSystem,
                      const interactor_ptr_t& interactor,
                      const sim_param_ptr_t& param)
    {
        // Setup
        static bool setup = false;
        static temperature_t temperature = param->get<real_t>("simulation.temperature");
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<real_t> disCoordinate(0.0, RANGE);
        static std::uniform_real_distribution<real_t> dis01(0.0, 1.0);
        static const real_t kT = units::mu<real_t>::KB * temperature();
        if ( !setup ) {      
            gen.seed(util::seedValue<std::size_t>());
            setup = true;
        }

        // Current position.
        position_t ri = particle->position();
        
        // Calculate initial (current) energy.
        auto result_i = interactor->interact(particle, particleSystem);
        auto energy_i = result_i.first + result_i.second;
        if ( energy_i() >= conf::LARGE || std::isnan(energy_i()) ) {
            energy_i = conf::LARGE;
        }

        // Displace/move particle.
        position_t rf;                         // Displaced position.
        for ( std::size_t k = 0; k != 3; ++k ) {
            rf[k] = ri[k] - 0.5 * RANGE + disCoordinate(gen);
        }
        particle->position(rf);
                
        // Calculate final (new) energy.
        auto result_f = interactor->interact(particle, particleSystem);
        auto energy_f = result_f.first + result_f.second;
        if ( energy_f()  >= conf::LARGE || std::isnan(energy_f()) ) {
            energy_f = 2.0 * conf::LARGE;
        }
        
        // Accept or reject displacement.
        energy_t difference = energy_f() - energy_i();
        real_t difference_over_kT = difference() / kT;
        if ( difference_over_kT < LIMIT && energy_f() < conf::LARGE) {
            if ( difference_over_kT > 0.0 ) {
                real_t w = std::exp(-difference_over_kT);
                real_t rv = dis01(gen);
                if ( rv > w ) {
                    // Reject. Restore previous position.
                    particle->position(ri);
                    return std::make_tuple(0.0, 0.0, false);
                } else {
                    // Accept. Keep new position.
                    return std::make_tuple(result_f.first - result_i.first,
                                           result_f.second - result_i.second,
                                           true);
                }
            } else {
                // Accept. Keep new position.
                return std::make_tuple(result_f.first - result_i.first,
                                       result_f.second - result_i.second,
                                       true);
            }
        } else {
            // Reject. Restore original position.
            particle->position(ri);
            return std::make_tuple(0.0, 0.0, false);
        }
    }

    static SimulationData 
    displaceOneParticle_(const std::vector<p_ptr_t>& all,
                         const p_system_ptr_t& particleSystem,
                         const interactor_ptr_t& interactor,
                         const sim_param_ptr_t& param)
    {
        static util::Logger logger("simploce::displaceOneParticle_");

        // Set up.
        static bool setup = false;
        temperature_t temperature = param->get<real_t>("simulation.temperature", 298.15);
        static energy_t bonded = 0.0;
        static energy_t nonBonded = 0.0;
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<std::size_t> dis(0, all.size() - 1);
        if ( !setup ) {
            gen.seed(util::seedValue<std::size_t>());
            setup = true;
        }

        // Randomly select one particle.
        auto index = dis(gen);
        auto particle = all[index];
        logger.debug("Selected particle with id '" + particle->id() + "'.");

        // Displace selected particle.
        auto result = displaceParticle_(particle, particleSystem, interactor, param);

        // Store some data.
        SimulationData data;
        bonded += std::get<0>(result);
        nonBonded += std::get<1>(result);
        data.bonded = bonded;
        data.nonBonded = nonBonded;
        data.accepted = std::get<2>(result);
        data.temperature = temperature;
        
        return std::move(data);
    }
    
    MonteCarlo::MonteCarlo(sim_param_ptr_t simulationParameters,
                           interactor_ptr_t interactor) :
        simulationParameters_{std::move(simulationParameters)}, interactor_{std::move(interactor)} {
    }

    SimulationData
    MonteCarlo::displace(const p_system_ptr_t& particleSystem) const {
        return particleSystem->doWithAll<SimulationData>([this, particleSystem] (
                const std::vector<p_ptr_t>& all) {
            return std::move(displaceOneParticle_(all,
                                                     particleSystem,
                                                     this->interactor_,
                                                     this->simulationParameters_));
        });
    }

}