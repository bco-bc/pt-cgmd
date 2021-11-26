/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 22 September 2019, 15:44
 */

#ifndef A_TYPES_HPP
#define A_TYPES_HPP

#include "simploce/util/param.hpp"
#include <memory>

namespace simploce {

    class Analyzer;
    class Gr;
    class DipoleMoment;

    /**
     * Analyzer pointer typeName.
     */
    using a_ptr_t = std::shared_ptr<Analyzer>;

    /**
     * Dipole moment analyzer pointer typeName.
     */
    using dm_ptr_t = std::shared_ptr<DipoleMoment>;

    /**
     * Gr pointer typeName.
     */
    using gr_ptr_t = std::shared_ptr<Gr>;

    /**
     * Analysis parameters.
     */
    using a_param_t = param::param_t;

    /**
     * Analysis parameters pointer typeName.
     */
    using a_param_ptr_t = std::shared_ptr<a_param_t>;
}

#endif /* A_TYPES_HPP */

