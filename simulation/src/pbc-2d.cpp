/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on February 1, 2024.
 */

#include <utility>

#include "simploce/simulation/pbc-2d.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/bc-util.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"

namespace simploce {

    PBC_2D::PBC_2D(simploce::box_ptr_t box,
                   Direction d1,
                   Direction d2,
                   Direction reinsert) :
        boundary_condition_impl{}, box_{std::move(box)}, directions_{d1.value(), d2.value()},
        reinsert_{reinsert.value()} {
    }

    dist_vect_t
    PBC_2D::apply(const position_t &ri,
                  const position_t &rj) const {
        static util::Logger logger{"simploce::PBC_2D::apply()"};
        logger.trace("Entering.");

        static int counter = 0;
        if (counter == 0) {
            logger.info("Applying periodicity in the x- andy-direction, but not in the z-direction.");
            counter += 1;
        }

        const box_t& box = *box_;
        auto rij = ri - rj;
        for (auto k : directions_) {
            auto box_k = box[k];
            auto dr = ri[k] - rj[k];
            auto ratio = dr / box_k;
            real_t n = util::nint(ratio);
            rij[k] = dr - n * box_k;
        }
        return rij;
    }

    velocity_t
    PBC_2D::apply(const velocity_t &v,
                  const position_t &r) const {
        return boundary_condition_impl::apply(v, r);
    }

    position_t
    PBC_2D::placeInside(const position_t &r_out) const {
        const box_t& box = *box_;
        position_t r_in{r_out};
        for (int k : directions_) {
            auto rk = r_in[k];
            auto box_k = box[k];
            auto ratio = rk / box_k;
            real_t n = util::nint(std::floor(ratio));
            r_in[k] -= n * box_k;
        }
        return r_in;
    }

    position_t
    PBC_2D::apply(const position_t &r) const {
        static auto box_k = box_->operator[](reinsert_);
        auto r_in = r;
        auto r_k = r_in[reinsert_];
        if ( bc::crossed(r_k, reinsert_, box_) ) {
            r_in[reinsert_] = util::randomUniform(0, box_k);
        }
        return r_in;
    }

    void
    PBC_2D::applyToVelocities(const simploce::pg_ptr_t &particleGroup) const {
        boundary_condition_impl::applyToVelocities(particleGroup);
    }
}
