/*
 * File:   protonatable-bead.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 11:37 AM
 */

#ifndef PROTONATABLE_HPP
#define PROTONATABLE_HPP

#include "simploce/types/u-types.hpp"
#include <cstdlib>
#include <map>

namespace simploce {

    /**
     * Capable of binding and releasing protons.
     */
    struct protonatable {

        ~protonatable() {}

        /**
         * "Protonates" this protonatable.
         */
        virtual void protonate() {}

        /**
         * "Deprotonates" this protonatable.
         */
        virtual void deprotonate() {}

        /**
         * Changes the protonation state of this entity.
         * @param values Values associated with protonation state.
         */
         virtual void protonate(const std::map<std::string, real_t>& values) {}

        /**
         * Inquires whether this entity is protonated.
         * @return Result.
         */
        virtual bool isProtonated() const = 0;

        /**
         * Returns charge.
         * @return Charge.
         */
        virtual charge_t charge() const = 0;

        /**
         * Returns mass.
         * @return Mass.
         */
        virtual mass_t mass() const = 0;

    };
}

#endif /* PROTONATABLE_HPP */

