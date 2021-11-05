/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 5, 2019, 2:24 PM
 */

#ifndef AT_DISPLACER_HPP
#define AT_DISPLACER_HPP

#include "displacer.hpp"
#include "sim-data.hpp"
#include "s-types.hpp"
#include "simploce/particle/atomistic.hpp"

namespace simploce {

    /**
     * Displaces atoms of an atomistic simulation model.
     */
    struct AtomisticDisplacer : public Displacer {
        
        virtual ~AtomisticDisplacer() {}
        
        /**
         * Displaces atoms.
         * @param param Simulation parameters.
         * @param at Atomistic particle model.
         * @return Simulation data (e.g. kinetic energy, temperature, etc).
         */
        virtual SimulationData displace(const sim_param_t& param, 
                                        const at_mod_ptr_t& at) const = 0;
        
    };
    
}


#endif /* AT_DISPLACER_HPP */

