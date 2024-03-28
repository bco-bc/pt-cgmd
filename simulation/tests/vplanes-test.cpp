//
// Created by juffer on 2/28/24.
//

#include "simploce/potentials/vplanes.hpp"
#include "simploce/potentials/vplane.hpp"
#include "simploce/potentials/uniform-surface-charge-density.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/simulation/pbc-2d.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/util.hpp"

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    auto box = factory::box(60.0, 60.0, 120.0);
    bc_ptr_t bc = factory::pbc_2d(box, Direction::X, Direction::Y);

    real_t eps_r = 78.5;
    VirtualPlanes vPlanes{box, bc, 1.0, eps_r};
    dist_t delta = 0.0;
    auto uniformSrfCgDensity = std::make_shared<UniformSurfaceChargeDensity>(-0.15, FlatSurface{}, eps_r, bc, delta, false);

    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto factory = simploce::factory::particleSystemFactory(catalog);
    std::string fn = "/wrk3/simulation/NaClNextToChargedSurface/no-external-field-sf-vplane/NaCl-cg-surface-1.ps";
    auto electrolyte = factory::particleSystem(fn, catalog, false);
    std::cout << electrolyte->numberOfParticles() << ": Number of ions." << std::endl;

    vPlanes.initialize(electrolyte);

    auto particles = electrolyte->doWithAll<std::vector<p_ptr_t>>([vPlanes] (const std::vector<p_ptr_t>& all) {
        return all;
    });

    auto Na = particles[0];
    std::cout << std::endl;
    std::cout << "Initial:" << std::endl;
    std::cout << "Interaction of one Na with one plane:" << std::endl;
    std::cout << util::to_string(VirtualPlanes::surfaceChargeDensity()) + ": Total surface charge density BEFORE." << std::endl;
    std::cout << "Particle BEFORE: " << *Na << std::endl;
    auto plane = VirtualPlanes::virtualPlanes()[10];
    std::cout << "Virtual plane: " << *plane << std::endl;
    auto result1 = plane->operator()(Na);
    std::cout << "Energy: " << result1.first << std::endl;
    std::cout << std::endl;

    auto r = Na->position();
    r[2] = plane->location()();
    Na->position(r);

    // Update.
    vPlanes.update(electrolyte);

    std::cout << std::endl;
    std::cout << "Final:" << std::endl;
    std::cout << "Interaction of one Na with one plane:" << std::endl;
    std::cout << util::to_string(VirtualPlanes::surfaceChargeDensity()) + ": Total surface charge density AFTER." << std::endl;
    std::cout << "Particle AFTER: " << *Na << std::endl;
    plane = VirtualPlanes::virtualPlanes()[10];
    std::cout << "Virtual plane: " << *plane << std::endl;
    auto result2 = plane->operator()(Na);
    std::cout << "Energy: " << result2.first << std::endl;
    std::cout << "Difference: " << result2.first - result1.first << std::endl;
    std::cout << std::endl;

    fn = "/wrk3/tests/vplanes.dat";
    std::ofstream stream;
    util::open_output_file(stream, fn);
    auto planes = VirtualPlanes::virtualPlanes();
    for (const auto& p: planes) {
        stream << p->location() << p->surfaceChargeDensity() << std::endl;
    }
    stream.close();
    std::clog << fn << ": Virtual planes data written to this file." << std::endl;

    /*
    fn = "/wrk3/tests/energy.dat";
    util::open_output_file(stream, fn);
    planes = VirtualPlanes::virtualPlanes();
    real_t dr = 0.1;
    auto n = size_t(box->lengthZ() / dr);
    for (auto k = 0; k != n; ++k) {
        auto z = k * dr;
        if (k == 0) z += 0.000001;
        r[2] = z;
        Na->position(r);
        energy_t energyNa{0};
        for (auto& p : planes) {
            auto result = p->operator()(Na);
            energyNa += result.first;
        }
        auto uniform = uniformSrfCgDensity->operator()(Na);
        auto total = energyNa + uniform.first;
        stream << r[2] << " " << energyNa << " " << uniform.first << " " << total << std::endl;
    }
    stream.close();
    std::cout << "Interaction energy Na with virtual planes, uniformSrfCgDensity surface charge density, and total energy "
                 "as function of z is written to '" << fn << "'." << std::endl;

    vPlanes.complete();
*/
}
