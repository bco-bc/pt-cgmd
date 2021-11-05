/*
 * Author: André H. Juffer.
 * Created on 27/10/2021, 14:39.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef FORCEFIELDS_LJ_HPP
#define FORCEFIELDS_LJ_HPP

#include "pair-potential.hpp"
#include "no-bc.hpp"
#include "simploce/types/u-types.hpp"
#include "simploce/util/map2.hpp"


namespace simploce {

    template<typename P, typename BC = no_bc>
    class LennardJones: public pair_potential<P> {
    public:

        /**
         * Type holding LJ interaction parameters, C12 and C6
         */
        using param_t = MatrixMap<std::string, std::pair<real_t, real_t>>;

        /**
         * Constructor
         * @param param Interaction parameters.
         */
        LennardJones(const param_t& param);

        /**
         * Particle pointer type.
         */
        using particle_ptr_t = typename pair_potential<P>::particle_ptr_t;

        energy_t energy(const particle_ptr_t& p1, const particle_ptr_t & p2) const override;

        energy_t force(particle_ptr_t& p1, particle_ptr_t& p2) const override;

    private:

        MatrixMap<std::string, std::pair<real_t, real_t>> param_;
        BC bc_;
    };

    template<typename P, typename BC>
    LennardJones<P,BC>::LennardJones(const param_t& param) : param_{param}, bc_{} {
    }

    template<typename P, typename BC>
    energy_t
    LennardJones<P,BC>::energy(const particle_ptr_t& p1, const particle_ptr_t & p2) const {
        auto name_i = p1->spec()->name();
        auto name_j = p2->spec()->name();
        auto key = std::make_pair(name_i, name_j);
        auto value = param_.get(key);
        auto r_i = p1->position();
        auto r_j = p2->position();
        dist_vect_t rij = bc_.apply(r_i, r_j);
        real_t Rij = norm<real_t>(rij);
        real_t Rij2 = Rij * Rij;
        real_t Rij6 = Rij2 * Rij2 * Rij2;
        real_t Rij12 = Rij6 * Rij6;
        real_t t1 = value.first / Rij12;    // C12
        real_t t2 = value.second / Rij6;    // C6
        real_t LJ = t1 - t2;                // kj/mol
        return std::move(energy_t{LJ});
    }

    template<typename P, typename BC>
    energy_t
    LennardJones<P,BC>::force(particle_ptr_t& p1, particle_ptr_t& p2) const {
        auto name_i = p1->spec()->name();
        auto name_j = p2->spec()->name();
        auto key = std::make_pair(name_i, name_j);
        auto value = param_.get(key);  // C12, C6

        auto r_i = p1->position();
        auto r_j = p2->position();
        dist_vect_t rij = bc_.apply(r_i, r_j);
        real_t Rij = norm<real_t>(rij);
        real_t Rij2 = Rij * Rij;
        real_t Rij6 = Rij2 * Rij2 * Rij2;
        real_t Rij12 = Rij6 * Rij6;
        real_t t1 = value.first / Rij12;    // C12
        real_t t2 = value.second / Rij6;    // C6
        real_t LJ = t1 - t2;                // Energy in kj/mol

        dist_vect_t uv = rij/Rij;                         // Unit vector
        real_t dLJdR = -6.0 * ( 2.0 * t1 - t2 ) / Rij;    // kJ/(mol nm)
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = -dLJdR * uv[k];    // kJ/(mol nm)
        }

        auto f1 = p1->force();
        f1 += f;
        p1->force(f1);
        auto f2 = p2->force();
        f2 -= f;
        p2->force(f2);

        return std::move(energy_t{LJ});

    }
}

#endif //FORCEFIELDS_LJ_HPP
