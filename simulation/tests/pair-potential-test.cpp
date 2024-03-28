/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#include "simploce/potentials/lj.hpp"
#include "simploce/potentials/halve-attractive-qp.hpp"
#include "simploce/potentials/halve-attractive-hp.hpp"
#include "simploce/potentials/hp.hpp"
#include "simploce/potentials/solid-sphere-dsf.hpp"
#include "simploce/potentials/gauss-sf.hpp"
#include "simploce/potentials/soft-repulsion.hpp"
#include "simploce/potentials/gauss-sf-sr.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/potentials/force-field.hpp"
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

void mc(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    box_ptr_t box = factory::box(conf::LARGE);
    bc_ptr_t bc = factory::pbc(box);
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

void testHAQP(const std::string& spec1,
              const std::string& spec2,
              const spec_catalog_ptr_t &catalog,
              const ff_ptr_t &forceField) {
    std::ofstream ostream;
    std::string fn = "/wrk3/tests/" + spec1 + "-" + spec2 + "-ha-qp.dat";
    util::open_output_file(ostream, fn);

    box_ptr_t box = factory::box(conf::LARGE);
    bc_ptr_t bc = factory::pbc(box);
    pair_potential_ptr_t ptr(new HalveAttractiveQP(forceField, bc));
    auto p1 = Bead::create("1", 0, "P1", catalog->lookup(spec1));
    p1->position(position_t{0.0, 0.0, 0.0});
    auto p2 = Bead::create("2", 1, "P2", catalog->lookup(spec2));
    p2->position(position_t{0.0, 0.0, 0.0});

    // Calculate potential.
    std::cout << "HalveAttractiveHQ:" << std::endl;
    std::cout << "Total " << spec1 << " - " << spec2 << " interaction energy and force on -second- particle." << std::endl;
    std::cout << "Results are written to '" << fn << "'." << std::endl;
    real_t dz = 0.01;
    dist_t rc = 2.6;

    p2->position(position_t{0.0, 0.0, 0.51});
    ptr->operator()(p1, p2);

    int n = util::nint(rc()/dz);
    for (int i = 1; i <= n; ++i) {
        real_t z = i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = ptr->operator()(p1, p2);
        // Force on second particle!
        auto f2 = -result.second[2];
        //std::cout << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
        ostream << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
    }
    ostream.close();
}

void testHAHP(const std::string& spec1,
              const std::string& spec2,
              const spec_catalog_ptr_t &catalog,
              const ff_ptr_t &forceField) {
    std::ofstream ostream;
    std::string fn = "/wrk3/tests/" + spec1 + "-" + spec2 + "-ha-hp.dat";
    util::open_output_file(ostream, fn);

    box_ptr_t box = factory::box(conf::LARGE);
    bc_ptr_t bc = factory::pbc(box);
    pair_potential_ptr_t ptr(new HalveAttractiveHP(forceField, bc));
    auto p1 = Bead::create("1", 0, "P1", catalog->lookup(spec1));
    p1->position(position_t{0.0, 0.0, 0.0});
    auto p2 = Bead::create("2", 1, "P2", catalog->lookup(spec2));
    p2->position(position_t{0.0, 0.0, 0.0});

    // Calculate potential.
    std::cout << "HalveAttractiveHP:" << std::endl;
    std::cout << "Total " << spec1 << " - " << spec2 << " interaction energy and force on -second- particle." << std::endl;
    std::cout << "Results are written to '" << fn << "'." << std::endl;
    real_t dz = 0.01;
    dist_t rc = 2.6;
    real_t z0 = 1.0e-10;
    int n = util::nint(rc()/dz);
    for (int i = 0; i <= n; ++i) {
        real_t z = z0 + i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = ptr->operator()(p1, p2);
        // Force on second particle!
        auto f2 = -result.second[2];
        //std::cout << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
        ostream << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
    }
    ostream.close();
}

void testHP(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    std::ofstream ostream;
    std::string fn = "/wrk3/tests/hp.dat";
    util::open_output_file(ostream, fn);
    auto box = factory::box(10, 10, 40);
    bc_ptr_t bc = factory::pbc(box);
    pair_potential_ptr_t ptr(new HP(forceField, bc));
    auto p1 = Bead::create("1", 0, "T1", catalog->lookup("PMU"));
    p1->position(position_t{0.0, 0.0, 0.0});
    auto p2 = Bead::create("1", 1, "T2", catalog->lookup("PMU"));
    p2->position(position_t{0.0, 0.0, 1.0});

    p2->position(position_t{0.0, 0.0, 0.51});
    ptr->operator()(p1, p2);

    dist_t rc = 40.0;
    real_t dz = 0.1;
    int n = util::nint(rc()/dz);
    for (int i = 1; i < n; ++i) {
        real_t z = i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = ptr->operator()(p1, p2);
        // Force on second particle!
        auto f2 = -result.second[2];
        std::cout << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
        ostream << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
    }
    ostream.close();
}

void testSolidSphere_DSF(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    std::ofstream ostream;
    std::string fn = "/wrk3/tests/solid-sphere-dsf.dat";
    util::open_output_file(ostream, fn);
    auto box = factory::box(1000, 1000, 1000);
    bc_ptr_t bc = factory::pbc(box);
    dist_t rc = 2.0;
    dist_t radius = 1.0;
    pair_potential_ptr_t ptr(new SolidSphere_DSF(rc, forceField, box, bc, radius, false));
    auto p1 = Bead::create("1", 0, "T1", catalog->lookup("CW"));
    p1->position(position_t{0.0, 0.0, 0.0});
    auto p2 = Bead::create("1", 1, "T2", catalog->lookup("CW"));
    p2->position(position_t{0.0, 0.0, 0.0});

    real_t dz = 0.1;
    int n = util::nint(rc()/dz);
    for (int i = 0; i < n+1; ++i) {
        real_t z = i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = ptr->operator()(p1, p2);
         // Force on second particle!
        auto f2 = -result.second[2];
        std::cout << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
        ostream << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
    }
}

void testGaussianSF(const std::string& spec1,
                    const std::string& spec2,
                    const spec_catalog_ptr_t &catalog,
                    const ff_ptr_t &forceField) {
    std::ofstream ostream;
    std::string fn = "/wrk3/tests/" + spec1 + "-" + spec2 + "-gaussian_sf.dat";
    util::open_output_file(ostream, fn);
    auto box = factory::box(1000, 1000, 1000);
    bc_ptr_t bc = factory::pbc(box);
    dist_t rcSR = 1.0;
    dist_t rcLR = 2.0 * rcSR;
    pair_potential_ptr_t ptr(new GaussianSF(forceField, bc, rcLR, true));
    auto p1 = Bead::create("1", 0, "P1", catalog->lookup(spec1));
    p1->position(position_t{0.0, 0.0, 0.0});
    auto p2 = Bead::create("1", 1, "P2", catalog->lookup(spec2));
    p2->position(position_t{0.0, 0.0, 0.0});

    std::cout << "GaussianSF:" << std::endl;
    std::cout << "Total " << spec1 << " - " << spec2 << " interaction energy and force on -second- particle." << std::endl;
    std::cout << "Results are written to '" << fn << "'." << std::endl;
    real_t dz = 0.01;
    int n = util::nint(2.0*rcSR() / dz);
    for (int i = 0; i < n+1; ++i) {
        real_t z = 0.0001 + i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto result = ptr->operator()(p1, p2);
         // Force on second particle!
        auto f2 = -result.second[2];
        //std::cout << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
        ostream << z << conf::SPACE << result.first << conf::SPACE << f2 << std::endl;
    }
}

void testSoftRepulsion(const std::string& spec1,
                       const std::string& spec2,
                       const spec_catalog_ptr_t &catalog,
                       const ff_ptr_t &forceField) {
    std::ofstream ostream;
    std::string fn = "/wrk3/tests/" + spec1 + "-" + spec2 + "-soft_repulsion.dat";
    util::open_output_file(ostream, fn);
    auto box = factory::box(1000, 1000, 1000);
    bc_ptr_t bc = factory::pbc(box);
    dist_t rcSR = 1.0;
    auto sr = std::make_shared<SoftRepulsion>(forceField, bc, rcSR);
    auto p1 = Bead::create("1", 0, "P1", catalog->lookup(spec1));
    p1->position(position_t{0.0, 0.0, 0.0});
    auto p2 = Bead::create("1", 1, "P2", catalog->lookup(spec2));
    p2->position(position_t{0.0, 0.0, 0.0});

    std::cout << "SoftRepulsion:" << std::endl;
    std::cout << "Total " << spec1 << " - " << spec2 << " interaction energy and force on -second- particle." << std::endl;
    std::cout << "Results are written to '" << fn << "'." << std::endl;
    real_t dz = 0.01;
    int n = util::nint(2.0 * rcSR() / dz);
    for (int i = 0; i < n+1; ++i) {
        real_t z = 0.0001 + i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto total = sr->operator()(p1, p2);
        // Force on second particle!
        auto f2 = -total.second[2];
        // std::cout << z << conf::SPACE << total.first << conf::SPACE << f2 << std::endl;
        ostream << z << conf::SPACE << total.first << conf::SPACE << f2 << std::endl;
    }

    ostream.flush();
    ostream.close();
}

void testGaussianSF_SoftRepulsion(const std::string& spec1,
                                  const std::string& spec2,
                                  const spec_catalog_ptr_t &catalog,
                                  const ff_ptr_t &forceField) {
    std::ofstream ostream;
    std::string fn = "/wrk3/tests/" + spec1 + "-" + spec2 + "-gaussian_sf_soft_repulsion.dat";
    util::open_output_file(ostream, fn);
    auto box = factory::box(1000, 1000, 1000);
    bc_ptr_t bc = factory::pbc(box);
    dist_t rcSR = 1.0;
    dist_t rcLR = 2.0 * rcSR;
    auto gaussSF = std::make_shared<GaussianSF>(forceField, bc, rcLR, true);
    auto sr = std::make_shared<SoftRepulsion>(forceField, bc, rcSR);
    auto ptr = std::make_shared<GaussianSF_SoftRepulsion>(forceField, bc, gaussSF, sr);
    auto p1 = Bead::create("1", 0, "P1", catalog->lookup(spec1));
    p1->position(position_t{0.0, 0.0, 0.0});
    auto p2 = Bead::create("1", 1, "P2", catalog->lookup(spec2));
    p2->position(position_t{0.0, 0.0, 0.0});

    std::cout << "GaussianSF+SoftRepulsion:" << std::endl;
    std::cout << "Total " << spec1 << " - " << spec2 << " interaction energy and force on -second- particle." << std::endl;
    std::cout << "Results are written to '" << fn << "'." << std::endl;
    real_t dz = 0.01;
    int n = util::nint(2.0 * rcSR() / dz);
    for (int i = 0; i < n+1; ++i) {
        real_t z = 0.0001 + i * dz;
        p2->position(position_t{0.0, 0.0, z});
        auto total = ptr->operator()(p1, p2);
        // Force on second particle!
        auto f2 = -total.second[2];
        // std::cout << z << conf::SPACE << total.first << conf::SPACE << f2 << std::endl;
        ostream <<  z << conf::SPACE << total.first << conf::SPACE << f2 << conf::SPACE << f2 << std::endl;
    }

    ostream.flush();
    ostream.close();
}

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGDEBUG);

    std::string fnSpecs = "/wrk3/tests/particle-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fnSpecs);
    auto catalog = factory::particleSpecCatalog(stream);
    stream.close();

    std::string fnInteractions = "/wrk3/tests/interaction-parameters.dat";
    // std::string fnInteractions = "/localdisk/resources/interaction-parameters-polymer-solution.dat";
    util::open_input_file(stream, fnInteractions);
    auto forceField = factory::forceField(stream, catalog);
    stream.close();

    //test1(catalog, forceField);

    //testHAQP("CW", "DP", catalog, forceField);
    //testHAHP("C-W", "D-P", catalog, forceField);

    // testSolidSphere_DSF(catalog, forceField);

    // testHP(catalog, forceField);

    //testSolidSphere_DSF(catalog, forceField);

    testGaussianSF("DP", "DP", catalog, forceField);
    testGaussianSF("CW", "DP", catalog, forceField);
    testGaussianSF_SoftRepulsion("CW", "CW", catalog, forceField);
    //testSoftRepulsion("C-W", "C-W", catalog, forceField);

    return EXIT_SUCCESS;
}