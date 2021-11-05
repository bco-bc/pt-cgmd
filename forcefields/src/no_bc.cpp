/*
 * Author: André H. Juffer.
 * Created on 27/10/2021, 15:46.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/forcefields/no-bc.hpp"

namespace simploce {

    dist_vect_t
    no_bc::apply(const position_t& r1, const position_t& r2) const {
        return std::move(r1 - r2);
    }
}