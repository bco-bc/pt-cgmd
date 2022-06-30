/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 14:24 PM
 */

#ifndef PARTICLES_CHEM_TYPES_HPP
#define PARTICLES_CHEM_TYPES_HPP

#include <memory>

namespace simploce {

    // Forward declarations.
    class content_handler;
    class chem_reader;
    class InputSource;

    /**
     * Content handler pointer type.
     */
    using cont_handler_ptr_t = std::shared_ptr<content_handler>;

    /**
     * Chemical reader pointer type.
     */
    using chem_reader_ptr_t = std::shared_ptr<chem_reader>;

    /**
     * Input source pointer type.
     */
    using input_source_ptr_t = std::shared_ptr<InputSource>;


}

#endif //PARTICLES_CHEM_TYPES_HPP
