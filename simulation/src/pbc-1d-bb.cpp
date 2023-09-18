/*
 * Author: André H. Juffer.
 * Created on 22/12/2021, 19:14.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/pbc-1d-bb.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/util/util.hpp"

namespace simploce {

    PBC_1D_BB::PBC_1D_BB(box_ptr_t box,
                         Direction direction) :
        box_{std::move(box)}, direction_{direction} {
    }

    dist_vect_t
    PBC_1D_BB::apply(const position_t& ri, const position_t& rj) const {
        static int index = direction_.value();
        const box_t& box = *box_;

        dist_vect_t rij = ri - rj;
        real_t box_k = box[index];
        real_t dr = ri[index] - rj[index];
        real_t ratio = dr / box_k;
        real_t n = util::nint(ratio);
        rij[index] = dr - n * box_k;

        return rij;
    }

    position_t
    PBC_1D_BB::placeInside(const position_t& r_out) const {
        PBC pbc(box_);
        return pbc.placeInside(r_out);
    }

    velocity_t
    PBC_1D_BB::apply(const velocity_t &v,
                     const position_t& r) const {

        // Here, index refers to the direction -no- PDB is applied.
        auto indices = [] (int index) {
            if (index == 0) {
                return std::vector<std::size_t>{1, 2};
            } else if (index == 1) {
                return std::vector<std::size_t>{0, 2};
            } else {
                return std::vector<std::size_t>{0, 1};
            }
        };

        static int index = direction_.value();
        static auto ks = indices(index);

        const box_t& box = *box_;
        bool crossed{false};
        for (auto k: ks) {
            auto box_k = box[k];
            if (r[k] > box_k || r[k] < 0.0)
                // Particle crossed the boundary in the k-direction.
                crossed = true;
        }
        return crossed ? -1.0 * v : v;
    }
}

