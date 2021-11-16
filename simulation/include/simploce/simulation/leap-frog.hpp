/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on August 16, 2019, 12:36 PM
 */

#ifndef LEAP_FROG_HPP
#define LEAP_FROG_HPP

#include "cg-displacer.hpp"
#include "at-displacer.hpp"
#include "sim-data.hpp"
#include "s-types.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"

namespace simploce {        
    
    /**
     * Standard leap frog algorithm for MD simulations.
     * @param M Particle model type.
     */
    class LeapFrog {
    public:    
        
        LeapFrog(sim_param_ptr_t simulationParameters);
                
        /**
         * Displaces atoms of an atomistic model.
         * @param at Atomistic model.
         * @return Kinetic energy, temperature.
         */
        SimulationData 
        displace(std::vector<std::shared_ptr<Particle>> &particles) const;
        
    private:

        sim_param_ptr_t simulationParameters_;
    };
    

}

#endif /* LEAP_FROG_HPP */

