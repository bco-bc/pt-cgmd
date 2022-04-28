/*
 * Author: André H. Juffer.
 * Created on 21/12/2021, 13:44.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/plane.hpp"
#include <utility>
#include <stdexcept>

namespace simploce {

    Plane::Plane(std::string value) :
        value_{std::move(value)} {
    }

    bool
    Plane::operator == (const Plane& plane) {
        return (value_ == plane.value_);
    }

    Plane
    Plane::valueOf(const std::string &value) {
        if (value == "yz") {
            return YZ;
        } else if (value == "zx") {
            return ZX;
        } else if (value == "xy") {
            return XY;
        } else {
            throw std::domain_error(value + ": No such plane.");
        }
    }

    Plane
    Plane::XY{"yx"};

    Plane
    Plane::YZ{"yz"};

    Plane
    Plane::ZX{"zx"};

    std::string
    Plane::toString() const {
        return "Plane(plane=" + value_ + ")";
    }

}

