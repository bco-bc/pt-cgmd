/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 14:24 PM
 */

#ifndef PARTICLES_CHEM_READER_HPPq
#define PARTICLES_CHEM_READER_HPP

#include "chem-types.hpp"

namespace simploce {

    /**
      * Interface for a reader that parses chemical content from an input source.
      * An implementation of this interface assumes a particular format of the chemical content.
    */
    struct chem_reader {

        virtual ~chem_reader();

        /**
         * Parses the content of the input source.
         * @param handler Content handler.
         * @param source Input source.
         */
        virtual void parse(const cont_handler_ptr_t& handler,
                           const input_source_ptr_t& source) = 0;

    };

}

#endif //PARTICLES_CHEM_READER_HPP
