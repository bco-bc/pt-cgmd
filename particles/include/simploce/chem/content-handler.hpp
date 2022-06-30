/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 15:01 PM
 */

#ifndef PARTICLES_CONTENT_HANDLER_HPP
#define PARTICLES_CONTENT_HANDLER_HPP

#include "simploce/types/u-types.hpp"
#include <string>

namespace simploce {

    /**
     * Interface for a content handler. The handler receives notification of chemical content.
     */
    struct content_handler {

        virtual ~content_handler();

        /**
         * Receive notification of beginning of chemical content
         * @param title Title of chemical content. Serves as an identification of the content.
         */
        virtual void start(const std::string &title) = 0;

        /**
         * Receive notification of end of chemical content
         */
        virtual void end() = 0;

        /**
         * Receive notification of the beginning of an atom group, e.g., residue.
         * @param name Name of atom group.
         */
        virtual void startAtomGroup(const std::string &name) = 0;

        /**
         * Receive notification of the end of an atom group.
         */
        virtual void endAtomGroup() = 0;

        /**
         * Receive notification of the beginning of a molecule
         * @param name Name of the molecule
         */
        virtual void startMolecule(const std::string &name) = 0;

        /**
         * Receive notification of the end of a molecule
         */
        virtual void endMolecule() = 0;

        /**
         * Receive notification of the beginning of an atom
         * @param name Name of the atom
         */
        virtual void startAtom(const std::string &name) = 0;

        /**
         * Receive notification of the atom coordinates.
         * @param r Position.
         */
        virtual void atomCoordinates(const position_t& r) = 0;

        /**
         * Receive notification of atom charge.
         * @param charge Charge value.
         */
        virtual void atomCharge(const charge_t& charge) = 0;

        /**
         * Receive notification of the end of an atom
         */
        virtual void endAtom() = 0;

        /**
         * Receives notification of an index (e.g. residue number, atom number).
         * @param index Index
         */
        virtual void index(int index) = 0;

    };

}

#endif //PARTICLES_CONTENT_HANDLER_HPP
