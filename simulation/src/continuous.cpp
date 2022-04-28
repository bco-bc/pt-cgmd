//
// Created by ajuffer on 10/28/21.
//

#include "simploce/simulation/continuous.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/units/units-mu.hpp"

namespace simploce {

    Continuous::Continuous(): x_{0.0}, dx_dt_{0.0}, protonated_{false} {
    }

    void
    Continuous::protonate() {
        // This should make dx_dt_ > 0.
        protonated_ = true;
    }

    void Continuous::deprotonate() {
         // This should make dx_dt_ < 0.
        protonated_ = false;
    }

    real_t
    Continuous::x() const {
        return x_;
    }

    void
    Continuous::x(real_t x) {
        x_ = x;
    }

    real_t
    Continuous::dx_dt() const {
    return dx_dt_;
    }

    void
    Continuous::dx_dt(real_t dx_dt) {
        dx_dt_ = dx_dt;
    }

    bool
    Continuous::isProtonated() const {
        return protonated_;
    }

    charge_t Continuous::charge() const {
        return x_ * units::mu<real_t>::PROTON_CHARGE;
    }

    mass_t
    Continuous::mass() const {
        return x_ * units::mu<real_t>::PROTON_MASS;
    }

    std::ostream& operator << (std::ostream& stream, const Continuous& continuous) {
        stream << conf::SPACE << continuous.x() << conf::SPACE << continuous.dx_dt();
        return stream;
    }

    std::istream& operator >> (std::istream& stream, Continuous& continuous) {
        real_t x, dx_dt;
        stream >> x >> dx_dt;
        continuous.x(x);
        continuous.dx_dt(dx_dt);
        return stream;
    }
}

