/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on November 22, 2021, 12:29 PM.
 */

#include "simploce/simulation/lj.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <string>

using namespace simploce;

int main() {
     util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    // Pair potential.
    std::string fileName = "/localdisk/resources/particles-specs.dat";
    auto catalog = factory::particleSpecCatalog(fileName);
    fileName = "/localdisk/resources/interaction-parameters.dat";
    auto forceField = factory::forceField(fileName, catalog);
    auto box = factory::box(2.0 * conf::CUTOFF_DISTANCE);
    auto bc = factory::boundaryCondition(box);
    LJ pairPotential{forceField, bc};

    // Particles.
    auto T1 = catalog->lookup("T1");
    p_ptr_t p1 = Atom::create("1", 0, "p1", T1);
    p1->position(position_t{0.0, 0.0, 0.0});
    p_ptr_t p2 = Atom::create("2", 1, "p2", T1);
    p2->position(position_t{0.0, 0.0, 0.0});

    // Calculate potential.
    real_t dz = 0.01;
    distance_t rc = properties::cutoffDistance(box);
    int n = util::nint(rc()/dz);
    for (int i = 1; i <= n; ++i) {
        real_t z = i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = pairPotential.operator()(p1, p2);
        // Force on second particle!
        auto f2 = -result.second[2];
        std::cout << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
    }

}