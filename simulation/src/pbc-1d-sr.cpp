/*
 * Author: André H. Juffer.
 * Created on 26/06/2023.
 *
 * Copyright (c) Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/pbc-1d-sr.hpp"
#include "simploce/simulation/bc-util.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/util/util.hpp"

namespace simploce {

    PBC_1D_SR::PBC_1D_SR(simploce::box_ptr_t box, simploce::Direction direction) :
        boundary_condition_impl{}, box_{std::move(box)}, direction_{direction} {
    }

    dist_vect_t
    PBC_1D_SR::apply(const simploce::position_t &ri,
                     const simploce::position_t &rj) const {
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
    PBC_1D_SR::placeInside(const simploce::position_t &r_out) const {
        static PBC pbc(box_);
        return pbc.placeInside(r_out);
    }

    velocity_t
    PBC_1D_SR::apply(const simploce::velocity_t &v,
                     const simploce::position_t &r) const {
        static auto nc = bc::normalComponents(direction_);

        auto vel = v;
        for (auto k: nc) {
            vel[k] = bc::crossed(r[k], k, box_) ? -vel[k] : vel[k];
        }
        return vel;
    }

    position_t
    PBC_1D_SR::apply(const simploce::position_t &r) const {
        return boundary_condition_impl::apply(r);
    }

    void
    PBC_1D_SR::applyToVelocities(const simploce::pg_ptr_t &particleGroup) const {
        static auto box = *box_;
        static auto nc = bc::normalComponents(direction_);

        auto r = particleGroup->position();
        std::vector<real_t> factor(3, 1.0);
        auto crossed = bc::crossed(r, box_, direction_);
        if (crossed) {
            for (auto k: nc) {
                factor[k] = bc::crossed(r[k], k, box_) ? -1.0 : 1.0;
            }
            for ( auto& p: particleGroup->particles() ) {
                auto v = p->velocity();
                for (auto k: nc) {
                    v[k] *= factor[k];
                }
                p->velocity(v);
            }
        }
    }

}
