/*
 * Author: André H. Juffer.
 * Created on 22/12/2021, 19:14.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/pbc-1d-bb.hpp"
#include "simploce/simulation/bc-util.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/util/util.hpp"

namespace simploce {

    PBC_1D_BB::PBC_1D_BB(box_ptr_t box,
                         Direction direction) :
        boundary_condition_impl{}, box_{std::move(box)}, direction_{direction} {
    }

    dist_vect_t
    PBC_1D_BB::apply(const position_t& ri,
                     const position_t& rj) const {
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
        static PBC pbc(box_);
        return pbc.placeInside(r_out);
    }

    position_t
    PBC_1D_BB::apply(const position_t &r) const {
        return boundary_condition_impl::apply(r);
    }

    velocity_t
    PBC_1D_BB::apply(const velocity_t &v,
                     const position_t& r) const {
        auto crossed = bc::crossed(r, box_, direction_);
        return crossed ? -1.0 * v : v;
    }

    void
    PBC_1D_BB::applyToVelocities(const simploce::pg_ptr_t &particleGroup) const {
        auto r = particleGroup->position();
        auto crossed = bc::crossed(r, box_, direction_);
        if (crossed) {
            for (auto& p: particleGroup->particles()) {
                auto v = p->velocity();
                v *= -1.0;
                p->velocity(v);
            }
        }
    }
}

