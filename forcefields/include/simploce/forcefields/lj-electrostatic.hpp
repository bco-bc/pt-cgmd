/*
 * Author: André H. Juffer.
 * Created on 27/10/2021, 15:00.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef FORCEFIELDS_LJ_ELECTROSTATIC_HPP
#define FORCEFIELDS_LJ_ELECTROSTATIC_HPP

#include "pair-potential.hpp"

namespace simploce {

    template <typename P>
    class LennardJonesElectrostatic : public pair_potential<P> {
    public:

        /**
         * Particle pointer type.
         */
        using particle_ptr_t = typename pair_potential<P>::particle_ptr_t;

        energy_t energy(const particle_ptr_t& p1, const particle_ptr_t & p2) const override;

        energy_t force(particle_ptr_t& p1, particle_ptr_t& p2) const override;

    };

    template<typename P>
    energy_t
    LennardJonesElectrostatic<P>::energy(const particle_ptr_t& p1, const particle_ptr_t & p2) const {
        return 0.0;
    }

    template<typename P>
    energy_t
    LennardJonesElectrostatic<P>::force(particle_ptr_t& p1, particle_ptr_t& p2) const {
        return 0.0;
    }
}

#endif //FORCEFIELDS_LJ_ELECTROSTATIC_HPP
