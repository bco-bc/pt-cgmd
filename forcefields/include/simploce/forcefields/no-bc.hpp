/*
 * Author: André H. Juffer.
 * Created on 27/10/2021, 15:39.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef FORCEFIELDS_NO_BC_HPP
#define FORCEFIELDS_NO_BC_HPP

#include "simploce/types/u-types.hpp"
#include <memory>

namespace simploce {

    /**
     * Applies no boundary condition.
     */
    struct no_bc {

        /**
         * Returns distance vector between r1 and r2.
         * @param r1 Position #1
         * @param r2 Position #2.
         * @return r1-r2.
         */
        dist_vect_t apply(const position_t& r1, const position_t& r2) const;

    };
}

#endif //FORCEFIELDS_NO_BC_HPP
