//
// Created by ajuffer on 2/1/24.
//

#include "simploce/potentials/uniform-surface-charge-density.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include <fstream>

using namespace simploce;

int main() {
    util::Logger logger("main");
    logger.changeLogLevel(util::Logger::LOGTRACE);

    auto box = factory::box(60.0, 60.0, 120.0);
    auto bc = factory::pbc_2d(box, Direction::X, Direction::Y);
    real_t eps = 78.5;
    dist_t delta = 0.0;
    auto potential = std::make_shared<UniformSurfaceChargeDensity>(-0.15, FlatSurface{}, eps, bc, delta, false);
    real_t dr = 0.1;
    int n = int(box->lengthZ() / dr);
    std::string fileName = "/wrk3/tests/uniform.dat";
    std::ofstream ostream;

    util::open_output_file(ostream, fileName);
    for (int i = 1; i != n; ++i) {
        auto z = i * dr;
        position_t r{0.0, 0.0, z};
        auto pair = potential->forceAndEnergy(r, 0.15, 1.0);
        std::cout << z << " " << pair.first << " " << pair.second << std::endl;
        ostream << z << " " << pair.first << " " << pair.second << std::endl;
    }
    ostream.close();
    std::cout << "Interaction energy and force written to \'" + fileName + "\'" << std::endl;
}
