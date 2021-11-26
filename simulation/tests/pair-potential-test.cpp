/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#include "simploce/simulation/lj.hpp"
#include "simploce/simulation/halve-attractive-qp.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/force-field.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/util.hpp"
#include <cstdlib>
#include <map>
#include <string>

using namespace simploce;

void test1(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    box_ptr_t box = factory::box(conf::LARGE);
    bc_ptr_t bc = factory::boundaryCondition(box);
    LJ lj{forceField, bc};
    pair_potential &pp = lj;

    std::map<std::string, pair_potential_ptr_t> map;
    auto spec = catalog->lookup("Ar");
    std::string key = spec->name() + "-" + spec->name();
    pair_potential_ptr_t ptr(new LJ(forceField, bc));
    auto pair = std::make_pair(key, ptr);
    map.emplace(pair);

    pair_potential &pairPotential = *map.at(key);
}

void testHAQP(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    box_ptr_t box = factory::box(conf::LARGE);
    bc_ptr_t bc = factory::boundaryCondition(box);
    pair_potential_ptr_t ptr(new HalveAttractiveQP(forceField, bc));
    auto p1 = Bead::create("1", 0, "CW", catalog->lookup("CW"));
    p1->position(position_t{0.0, 0.0, 0.0});
    auto p2 = Bead::create("1", 1, "DP", catalog->lookup("DP"));
    p2->position(position_t{0.0, 0.0, 0.0});

    // Calculate potential.
    // Calculate potential.
    real_t dz = 0.01;
    distance_t rc = properties::cutoffDistance(box);
    std::cout << "Cutoff distance: " << rc << std::endl;

    p2->position(position_t{0.0, 0.0, 0.51});
    ptr->operator()(p1, p2);

    int n = util::nint(rc()/dz);
    for (int i = 1; i <= n; ++i) {
        real_t z = i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = ptr->operator()(p1, p2);
        // Force on second particle!
        auto f2 = -result.second[2];
        std::cout << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
    }


}

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fnSpecs = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fnSpecs);
    auto catalog = factory::particleSpecCatalog(stream);
    stream.close();

    std::string fnInteractions = "/localdisk/resources/interaction-parameters.dat";
    util::open_input_file(stream, fnInteractions);
    auto forceField = factory::forceField(stream, catalog);
    stream.close();

    test1(catalog, forceField);

    testHAQP(catalog, forceField);

    return EXIT_SUCCESS;
}