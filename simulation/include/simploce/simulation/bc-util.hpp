/*
 * Author: André H. Juffer.
 * Created on 26/06/2023.
 *
 * Copyright (c) Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_BC_UTIL_HPP
#define SIMULATION_BC_UTIL_HPP

#include "../types/s-types.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {
    namespace bc {

        /**
         * Returns the other two normal directions parallel to the given direction that are used to check whether
         * a particle crossed the box boundary in these directions.
         * @param direction Direction.
         * @return Indices normal components.
         */
        std::vector<size_t>
        normalComponents(const Direction& direction);

        /**
         * Checks whether a particle at position r crossed the box boundary in the directions other than
         * the given direction.
         * @param r Position.
         * @param box Box
         * @param direction Direction.
         * @return True if crossing occurred, otherwise false.
         */
        bool
        crossed(const position_t &r, const box_ptr_t &box, const Direction &direction);

        /**
         * Checks whether a position vector component crossed the box boundary in the direction specified by index.
         * @param rc Position vector component.
         * @param index Index, specifies direction.
         * @param box Box.
         * @return True if crossing occurred, otherwise false.
         */
        bool
        crossed(real_t rc, size_t index, const box_ptr_t &box);

    }
}

#endif //SIMULATION_BC_UTIL_HPP
