/*
 * File:   displacer.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 13, 2019, 3:30 PM
 */

#ifndef DISPLACER_HPP
#define DISPLACER_HPP

#include <string>

namespace simploce {
    
    /**
     * Displaces or changes the state of simulation model.
     */
    struct Displacer {
        
        virtual ~Displacer() {}
        
        /**
         * Returns identifying name.
         * @return Identifying name.
         */
        virtual std::string id() const = 0;
    };
}

#endif /* DISPLACER_HPP */

