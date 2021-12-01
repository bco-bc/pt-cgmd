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

    class analyzer;
    class Gr;
    class DipoleMoment;
    class Distance;

    /**
     * Analyzer pointer type.
     */
    using a_ptr_t = std::shared_ptr<analyzer>;

    /**
     * Dipole moment analyzer pointer type.
     */
    using dm_ptr_t = std::shared_ptr<DipoleMoment>;

    /**
     * g(r) analyzer pointer type.
     */
    using gr_ptr_t = std::shared_ptr<Gr>;

    /**
     * Distance analyzer pointer type.
     */
    using d_ptr_t = std::shared_ptr<Distance>;

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

