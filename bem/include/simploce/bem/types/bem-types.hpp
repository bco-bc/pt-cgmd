/*
 * Author: André H. Juffer.
 * Created on 23/05/2022, 21:26.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef BEM_BEM_TYPES_HPP
#define BEM_BEM_TYPES_HPP

#include "simploce/types/u-types.hpp"
#include "simploce/surface/types/srf-types.hpp"
#include "simploce/util/param.hpp"
#include <memory>

namespace simploce {

    // Forward declarations.
    class bem_calculator;

    /**
     * Surface pointer type.
     */
    using surface_ptr_t = polyhedron_ptr_t;

    /**
     * Parameters type.
     */
    using param_t = param::param_t;
    using param_ptr_t = param::param_ptr_t;

    /**
     * Pointer type for the BEM calculator.
     */
    using bem_calc_ptr_t = std::shared_ptr<bem_calculator>;
}

#endif //BEM_BEM_TYPES_HPP
