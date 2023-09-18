/*
 * Author: André H. Juffer.
 * Created on 21/12/2021, 15:00.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/direction.hpp"
#include <stdexcept>

namespace simploce {

    Direction::Direction(char value) :
        value_(value) {
    }

    bool Direction::operator == (const Direction& direction) const {
        return value_ == direction.value_;
    }

    int Direction::value() const {
        if (value_ == 'x') {
            return 0;
        } else if (value_ == 'y') {
            return 1;
        } else {
            return 2;
        }
    }

    Direction
    Direction::X{'x'};

    Direction
    Direction::Y{'y'};

    Direction
    Direction::Z{'z'};

    Direction
    Direction::valueOf(char value) {
        if (value == 'x') {
            return Direction::X;
        } else if (value == 'y') {
            return Direction::Y;
        } else if (value == 'z') {
            return Direction::Z;
        } else {
            throw std::domain_error(value + ": no such direction.");
        }
    }
}

