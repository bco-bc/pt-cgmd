/*
 * Author: André H. Juffer.
 * Created on 21/12/2021, 13:55.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/flat-surface.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/util/util.hpp"
#include <cmath>

namespace simploce {

    FlatSurface::FlatSurface(Plane plane, const dist_t &distanceToPlane) :
            plane_{plane}, distanceToPLane_{distanceToPlane}, unitVector_{} {
        if (plane_ == Plane::XY) {
            unitVector_ = dist_vect_t{0.0, 0.0, 1.0};
            coordinate_ = 2;
        } else if (plane == Plane::YZ) {
            unitVector_ = dist_vect_t{1.0, 0.0, 0.0};
            coordinate_ = 0;
        } else {
            unitVector_ = dist_vect_t{0.0, 1.0, 0.0};
            coordinate_ = 1;
        }
    }

    std::pair<dist_t, dist_vect_t>
    FlatSurface::distanceTo(const position_t &r) const {
        dist_t R{std::fabs(distanceToPLane_() - r[coordinate_])};
        dist_vect_t ris = R() > distanceToPLane_() ? R() * unitVector_ : -R() * unitVector_;
        return std::make_pair(R, ris);
    }

    dist_vect_t
    FlatSurface::unitVectorPerpendicularTo() const {
        return unitVector_;
    }

    std::string
    FlatSurface::toString() const {
        return "FlatSurface(" + plane_.toString() +
        ", distanceToPlane=" + util::toString(distanceToPLane_()) +
        ", unitVector=" + util::toString(unitVector_) +
        ", coordinate=" + util::toString(coordinate_) +
        ")";
    }

    std::ostream& operator <<(std::ostream& stream, const FlatSurface& flatSurface) {
        stream << flatSurface.plane_.toString() << conf::SPACE << flatSurface.distanceToPLane_();
        stream << conf::SPACE << flatSurface.unitVector_ << conf::SPACE << flatSurface.coordinate_;
        return stream;
    }

}