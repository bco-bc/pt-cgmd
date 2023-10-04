/*
 * Author: André H. Juffer.
 * Created on 26/06/2023.
 *
 * Copyright (c) Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/bc-util.hpp"

namespace simploce {
    namespace bc {

        std::vector<size_t>
        normalComponents(const Direction& direction) {
            auto index = direction.value();
            if (index == 0) {
                return std::vector<size_t>{1, 2};
            } else if (index == 1) {
                return std::vector<std::size_t>{0, 2};
            } else {
                return std::vector<std::size_t>{0, 1};
            }
        }

        bool
        crossed(const position_t &r, const box_ptr_t &box, const Direction &direction) {
            static box_t& b = *box;
            static auto ks = normalComponents(direction);

            bool crossed{false};
            for (auto k: ks) {
               if ( bc::crossed(r[k], k, box) )
                    // Particle crossed the boundary in the k-direction.
                    crossed = true;
            }
            return crossed;
        }

        inline bool
        crossed(real_t rc, size_t index, const box_ptr_t &box) {
            static box_t& b = *box;
            auto box_k = b[index];
            return rc > box_k || rc < 0.0;
        }
    }
}
