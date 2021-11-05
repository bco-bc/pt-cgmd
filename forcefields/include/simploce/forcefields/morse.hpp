/*
 * Author: André H. Juffer.
 * Created on 25/10/2021, 21:45.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef FORCEFIELDS_MORSE_HPP
#define FORCEFIELDS_MORSE_HPP

#include "dw.hpp"
#include <cmath>

namespace simploce {

    /**
     * Double well potential U(x) constructed from two Morse potentials, where x is the location
     * in the double well potential.
     * @tparam P Particle type.
     */
    template <typename P>
    class Morse: public dw<P> {
    public:

        /**
         * Constructor. The default values result in barrier height of about 1 kT (2 kT) for the
         * particle moving from well #1 (#2) to well #2 (#1)
         * @param r0 Controls the location of the minima of the double well. The double well
         * possesses two minima (#1 and #2), one at r0 and another one at
         * r0 + R, respectively.
         * @param R Distance (in nm) between the two well locations.
         * @param wellDepth_1 Depth of well #1 in kT.
         * @param wellDepth_2 Depth of well #2 in kT.
         * @param T Temperature in K.
         */
        Morse(const distance_t& r0 = 0.1, const distance_t& R = 0.27,
                       real_t wellDepth_1 = 1.3, real_t fc_1 = 5000.0,
                       real_t wellDepth_2 = 2.3, real_t fc_2 = 5000.0,
                       const temperature_t& T = 298.0);

        /**
         * Particle pointer type.
         */
        using particle_ptr_t = typename dw<P>::particle_ptr_t;

        real_t U(real_t x) const override;

        real_t dUdx(real_t x) const override;

        real_t dU2dx2(real_t x) const override;

        /**
         * Energy.
         * @param p Particle located on the x-axis.
         * @return Energy.
         */
        energy_t energy(const particle_ptr_t& p) const override;

        energy_t force(particle_ptr_t& p) const override;

        void changeR(const distance_t& R);

    private:

        position_t r1_;
        distance_t r0_;
        distance_t R_;
        real_t wellDepth_1_;
        real_t fc_1_;
        real_t a_1_;
        real_t wellDepth_2_;
        real_t fc_2_;
        real_t a_2_;
        temperature_t T_;
    };

    template <typename P>
    Morse<P>:: Morse(const distance_t& r0, const distance_t& R,
                       real_t wellDepth_1, real_t fc_1,
                       real_t wellDepth_2, real_t fc_2,
                       const temperature_t& T) :
        dw<P>{}, r0_{r0}, R_{R},
        wellDepth_1_{wellDepth_1}, fc_1_{fc_1}, wellDepth_2_{wellDepth_2}, fc_2_{fc_2}, T_{T} {
        a_1_ = std::sqrt(fc_1_ / (2.0 * wellDepth_1_));
        a_2_ = std::sqrt(fc_2_ / (2.0 * wellDepth_2_));
    }

    template <typename P>
    real_t Morse<P>::U(real_t x) const {
        return 0.0;
    }

    template <typename P>
    real_t Morse<P>::dUdx(real_t x) const {
        return 0.0;
    }

    template <typename P>
    real_t Morse<P>::dU2dx2(real_t x) const {
        return 0.0;
    }

    template <typename P>
    energy_t Morse<P>::energy(const particle_ptr_t& p) const {
        real_t x = p->position()[0];
        return energy_t{U(x)};
    }

    template <typename P>
    energy_t
    Morse<P>::force(particle_ptr_t& p) const {
        auto energy = this->energy(p);
        return energy;
    }

    template <typename P>
    inline void
    Morse<P>::changeR(const distance_t &R) {
        R_ = R;
    }

}

#endif //FORCEFIELDS_MORSE_HPP
