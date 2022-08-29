/*
 * Author: André H. Juffer.
 * Created on 12/07/2022.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef BEM_RHS_HPP
#define BEM_RHS_HPP

#include "types/bem-types.hpp"
#include "bem-data.hpp"
#include <vector>

namespace simploce {
    namespace rhs {

        /**
         * Returns the right hand side of
         * @param nodes Nodes.
         * @param rcg Charge positions.
         * @param cg Charge values.
         * @param epsSolvent Dielectric constant solvent.
         * @param epsSolute Dielectric constant solute.
         * @return
         */
        BEMData::vector_t rightHandSide(std::vector<BEMData::node_t>& nodes,
                                        std::vector<position_t> &rcg,
                                        std::vector<charge_t> &cg,
                                        real_t epsSolvent,
                                        real_t epsSolute);
    }
}

#endif //BEM_RHS_HPP
