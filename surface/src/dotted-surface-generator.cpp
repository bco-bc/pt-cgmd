/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#include "simploce/surface/dotted-surface-generator.hpp"
#include "simploce/surface/nsc.hpp"
#include "simploce/surface/conf/srf-conf.hpp"
#include "simploce/util/math-constants.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <cmath>

namespace simploce {

    static std::pair<std::vector<dot_t>, area_t>
    nsc_(const std::vector<position_t>& positions, const std::vector<radius_t>& radii) {
        util::Logger logger("simploce::nsc_()");
        logger.trace("Entering.");

        // Convert to format expected by NSC.
        auto N = positions.size();
        auto *r = new double[3 * N];
        auto *rad = new double[N];
        for (std::size_t i = 0; i != N; ++i) {
            auto p = positions[i];
            rad[i] = radii[i]() + conf::H2O_RADIUS;
            auto i3 = 3 * i;
            for (std::size_t k = 0; k != 3; ++k) {
                r[i3 + k] = p[k];
            }
        }

        // Compute surface points.
        double *areas{nullptr};
        double *surfacePoints{nullptr};
        double area{0.0}, volume{0.0};
        int numberOfSurfacePoints{0};
        int result =
            NSC(r, rad, int(N), conf::DOT_DENSITY, FLAG_DOTS, &area, &areas, &volume, &surfacePoints, &numberOfSurfacePoints);
        logger.debug(std::to_string(result) + ": NSC return value.");
        logger.debug(std::to_string(area) + ": Dotted surface area.");
        logger.debug(std::to_string(volume) + ": Dotted surface enclosed volume.");
        logger.debug(std::to_string(numberOfSurfacePoints) + ": Number of generated surface points (dots).");

        // Convert back.
        std::vector<dot_t> dots;
        for (std::size_t i = 0; i != numberOfSurfacePoints; ++i) {
            auto i3 = 3 * i;
            real_t x = surfacePoints[i3];
            real_t y = surfacePoints[i3 + 1];
            real_t z = surfacePoints[i3 + 2];
            dots.emplace_back(dot_t{x, y, z});
        }

        // Free memory.
        delete[](r);
        delete[](rad);
        delete[](areas);
        delete(surfacePoints);

        // Done.
        logger.trace("Leaving.");
        return std::move(std::make_pair(dots, area));
    }

    std::pair<std::vector<dot_t>, area_t>
    dotted_surface_generator::spherical(radius_t radius, real_t density) {
        std::vector<dot_t> dots;
        auto area = 4.0 * math::constants<real_t>::PI * radius() * radius();
        auto numberOfDots = std::size_t(density * area);
        auto r = radius();
        auto r2 = 2.0 * r;
        for (std::size_t k = 0; k != numberOfDots; ++k) {
            auto v = util::randomUniform();
            auto x = -r + v * r2;

            v = util::randomUniform();
            auto y = -r + v* r2;

            v = util::randomUniform();
            auto z = -r + v * r2;

            dot_t dot{x, y, z};
            auto length = norm<real_t>(dot);
            dot *= (radius() / length);

            dots.emplace_back(dot);
        }
        return std::move(std::make_pair(dots, area));
    }

    std::pair<std::vector<dot_t>, area_t>
    dotted_surface_generator::general(std::vector<position_t>& positions,
                                      std::vector<radius_t>& radii) {
        return std::move(nsc_(positions, radii));
    }
}
