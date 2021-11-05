//
// Created by ajuffer on 10/28/21.
//

#include "simploce/simulation/discrete.hpp"
#include "simploce/units/units-mu.hpp"

namespace simploce {

    Discrete::Discrete() :  x_{0} {
    }

    void
    Discrete::protonate() {
        x_ = 1;
    }

    void
    Discrete::deprotonate() {
        x_ = 0;
    }

    bool
    Discrete::isProtonated() const {
        return x_ == 1;
    }

    charge_t Discrete::charge() const {
        return x_ * units::mu<real_t>::PROTON_CHARGE;
    }

    mass_t
    Discrete::mass() const {
        return x_ * units::mu<real_t>::PROTON_MASS;
    }

}

