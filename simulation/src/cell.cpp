/*
 * Author: André H. Juffer.
 * Created on 31/08/22, 20:35.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/cell.hpp"

namespace simploce {

    std::string
    Cell::stringLocation(const simploce::Cell::location_t &location) {
        auto i = std::get<0>(location);
        auto j = std::get<1>(location);
        auto k = std::get<2>(location);
        return std::to_string(i) + "-" + std::to_string(j) + "-" + std::to_string(k);
    }

    Cell::Cell(Cell::location_t location, position_t r) :
            particles_{}, location_{std::move(location)}, r_{r},
            stringLocation_{Cell::stringLocation(location_)} {
    }

    bool Cell::operator == (const Cell& cell) const {
        return location_ == cell.location_;
    }

    void
    Cell::assign(const simploce::p_ptr_t &particle) {
        particles_.emplace_back(particle);
    }

    const std::vector<p_ptr_t>&
    Cell::particles() const {
        return particles_;
    }

    Cell::location_t
    Cell::location() const {
        return location_;
    }

    std::string
    Cell::stringLocation() const {
        return stringLocation_;
    }

    position_t
    Cell::position() const {
        return r_;
    }

    void
    Cell::clear_() {
        particles_.clear();
    }


}
