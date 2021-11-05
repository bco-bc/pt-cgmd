/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 13, 2019, 3:22 PM
 */

#ifndef FORCEFIELD_HPP
#define FORCEFIELD_HPP

#include "s-types.hpp"
#include "simploce/util/map2.hpp"
#include <string>
#include <utility>

namespace simploce {
    
    /**
     * Force field.
     */
    struct ForceField {
        
        /**
         * Type for holding pairs (C12, C6) of LJ interaction parameters 
         * for pairs of particle specs, each identified by a specification name. 
         */
        using lj_params_t = MatrixMap<std::string, std::pair<real_t, real_t>>;
        
        /**
         * Type for holding parameters for electrostatic interactions.
         */
        using el_params_t = std::map<std::string, real_t>;
        
        virtual ~ForceField() {}
        
        /**
         * Returns identifying name.
         * @return Identifying name.
         */
        virtual std::string id() const = 0;
        
        /**
         * Returns Lennard-Jones and electrostatic interaction parameters.
         * @return Parameters,
         */
        virtual std::pair<lj_params_t, el_params_t> parameters() const = 0;
    };
}

#endif /* FORCEFIELD_HPP */

