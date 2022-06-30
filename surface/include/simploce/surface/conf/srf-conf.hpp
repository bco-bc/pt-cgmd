/*
 * Author: André H. Juffer.
 * Created on 19/05/2022, 20:32.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SURFACE_SRF_CONF_HPP
#define SURFACE_SRF_CONF_HPP

#include "simploce/types/u-types.hpp"

namespace simploce {
    namespace conf {

        /**
         * NSC surface dot density.
         */
        const int DOT_DENSITY = 1000.0;  // Per spherical object.

        const real_t H2O_RADIUS = 0.14;  // Radius of a water molecule (nm).

    }

}

#endif //SURFACE_SRF_CONF_HPP
