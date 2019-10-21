/*
 * The MIT License
 *
 * Copyright 2019 juffer.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "simploce/simulation/mc.hpp"
#include "simploce/simulation/sim-model.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/mu-units.hpp"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <memory>
#include <utility>
#include <cmath>

namespace simploce {
    
    static const real_t RANGE = 0.1;
    static const real_t LIMIT = 75.0;
    static const real_t LARGE = 1.0e+20;
    
    template <typename P>
    static std::pair<energy_t, bool> 
    displaceParticle_(std::shared_ptr<P>& particle,
                      const cg_sim_model_ptr_t& sm,
                      const sim_param_t& param)
    {
        // Setup
        static bool setup = false;
        static temperature_t temperature = param.get<real_t>("temperature");
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<real_t> disCoordinate(0.0, RANGE);
        static std::uniform_real_distribution<real_t> dis01(0.0, 1.0);
        static const real_t kT = MUUnits<real_t>::KB * temperature();
        if ( !setup ) {      
            gen.seed(util::seedValue<std::size_t>());
            setup = true;
        }

        // Current position.
        position_t ri = particle->position();
        
        // Calculate initial (current) energy.
        energy_t energy_i = sm->interact(particle, param).epot;
        if ( energy_i() >= LARGE || std::isnan(energy_i()) ) {
            energy_i = LARGE;
        }

        // Move particle.
        position_t rf;                         // Displaced position.
        for ( std::size_t k = 0; k != 3; ++k ) {
            rf[k] = ri[k] - 0.5 * RANGE + disCoordinate(gen);
        }
        particle->position(rf);
                
        // Calculate final (new) energy.
        energy_t energy_f = sm->interact(particle, param).epot;
        if ( energy_f()  >= LARGE || std::isnan(energy_f()) ) {
            energy_f = 2.0 * LARGE;
        }
        
        // Accept or reject displacement.
        energy_t difference = energy_f() - energy_i();
        real_t difference_over_kT = difference() / kT;
        if ( difference_over_kT < LIMIT && energy_f() < LARGE) {
            if ( difference_over_kT > 0.0 ) {
                real_t w = std::exp(-difference_over_kT);
                real_t rv = dis01(gen);
                if ( rv > w ) {
                    // Reject. Restore previous position.
                    particle->position(ri);
                    return std::make_pair(energy_t(), false);
                } else {
                    // Accept. Keep new position.
                    return std::make_pair(difference, true);
                }
            } else {
                // Accept. Keep new position.
                return std::make_pair(difference, true);;
            }
        } else {
            // Reject. Restore position.
            particle->position(ri);
            return std::make_pair(energy_t(), false);;
        }
    }
    
    template <typename P>
    static SimulationData 
    displaceOneParticle_(const std::vector<std::shared_ptr<P>>& all,
                         const cg_sim_model_ptr_t& sm,
                         const sim_param_t& param)
    {
        
        // Set up.
        static bool setup = false;        
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<std::size_t> dis(0, all.size() - 1);
        if ( !setup ) {
            gen.seed(util::seedValue<std::size_t>());
            setup = true;
        }
        
        auto index = dis(gen);
        auto particle = all[index];
        auto result = displaceParticle_(particle, sm, param);
        
        SimulationData data{};
        data.epot = result.second;
        data.accepted = result.second;
        
        return data;
    }
    
    MC<Bead>::MC(const cg_sim_model_ptr_t& sm) : sm_{sm}
    {        
    }
    
    void 
    MC<Bead>::perform(const sim_param_t& param,
                      std::ofstream& trajStream,
                      std::ofstream& dataStream)
    {
        const auto width = conf::WIDTH;
        const auto space = conf::SPACE;
        
        if ( sm_->size() == 0 ) {
            throw std::domain_error(
                "No particles! Nothing to simulate."
            );
        }

        std::size_t nsteps = param.get<std::size_t>("nsteps", 10000);
        std::size_t nwrite = param.get<std::size_t>("nwrite", 10);
        
        for (std::size_t counter = 1; counter <= nsteps; ++counter) {
            SimulationData data = 
                sm_->doWithAll<SimulationData>([this, param] (std::vector<bead_ptr_t>& all) {
                    return displaceOneParticle_<Bead>(all, this->sm_, param);
                });
            if ( counter % nwrite == 0 ) {
                dataStream << std::setw(width) << counter << space << data << std::endl;
                sm_->saveState(trajStream);
                trajStream.flush();
                dataStream.flush();
            }
        }
        
    }
    
}