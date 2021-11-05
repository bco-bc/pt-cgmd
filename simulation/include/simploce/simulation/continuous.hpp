//
// Created by ajuffer on 10/28/21.
//

#ifndef SIMULATION_CONTINUOUS_HPP
#define SIMULATION_CONTINUOUS_HPP

#include "simploce/particle/protonatable.hpp"
#include <iostream>

namespace simploce {

    /**
     * Represents a continuously changing protonation state. Its value is in [0,1].
     */
    class Continuous: public protonatable {
    public:

        Continuous();

        void protonate() override;

        void deprotonate() override;

        /**
         * Returns value of x.
         * @return x.
         */
        real_t x() const;

        /**
         * Sets x.
         * @param x Value, must be in interval [0,1].
         */
        void x(real_t x);

        /**
         * Returns value of dx/dt.
         * @return dx/dt.
         */
        real_t dx_dt() const;

        /**
         * Sets dx/dt, the "current".
         * @param dx_dt Value.
         */
        void dx_dt(real_t dx_dt);

        bool isProtonated() const override;

        charge_t charge() const override;

        mass_t mass() const override;

    private:

        real_t x_;
        real_t dx_dt_;
        bool protonated_;
    };

    std::ostream& operator << (std::ostream& stream, const Continuous& continuous);

    std::istream& operator >> (std::istream& stream, Continuous& continuous);
}

#endif //SIMULATION_CONTINUOUS_HPP
