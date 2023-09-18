/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on Jun 2, 2022
 */

#include "simploce/potentials/hp-sr.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/file.hpp"
#include <string>
#include <limits>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fileName = "/localdisk/resources/particles-specs.dat";
    auto catalog = factory::particleSpecCatalog(fileName);

    fileName = "/localdisk/resources/interaction-parameters-polymer-solution.dat";
    auto forceField = factory::forceField(fileName, catalog);

    std::clog << "Force field:" << std::endl;
    std::clog << *forceField << std::endl;

    dist_t cutoff = 1.0;
    auto box = factory::box(10.0, 10.0, 40.0);
    auto bc = factory::pbc(box);
    auto harmonic = std::make_shared<HP>(forceField, bc);
    auto softRepulsion = std::make_shared<SoftRepulsion>(forceField, bc, cutoff);
    HarmonicSoftRepulsion pairPotential(forceField, bc, softRepulsion, harmonic);

    // Particles. Should form a bond.
    // These specs should exist in bonds.
    auto spec_1 = catalog->lookup("PMU");
    auto spec_2 = catalog->lookup("PMU");
    p_ptr_t p1 = Bead::create("1", 0, "p1", spec_1);
    p1->position(position_t{0.0, 0.0, 0.0});
    p_ptr_t p2 = Bead::create("2", 1, "p2", spec_2);
    auto z0 = real_t(std::numeric_limits<float>::min());
    p2->position(position_t{0.0, 0.0, z0});

    // Calculate potential.
    fileName = "/wrk3/tests/hp-sr-potential.dat";
    std::ofstream ostream;
    util::open_output_file(ostream, fileName);
    real_t dz = 0.01;
    int n = util::nint(cutoff() / dz);
    n *= 4;
    for (int i = 1; i <= n; ++i) {
        real_t z = z0 + i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = pairPotential.operator()(p1, p2);
        std::cout << z << conf::SPACE << result.first << conf::SPACE << result.second << std::endl;
        ostream << z << conf::SPACE << result.first << conf::SPACE << result.second << std::endl;
    }
    ostream.close();
}
