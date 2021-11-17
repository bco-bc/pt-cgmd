/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#include "simploce/simulation/lj.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include <cstdlib>
#include <map>
#include <string>

using namespace simploce;

void test1(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    box_ptr_t box = factory::box(1.0);
    bc_ptr_t bc = factory::pbc(box);
    LJ lj{forceField, box, bc};
    pair_potential &pp = lj;

    std::map<std::string, pair_potential_ptr_t> map;
    auto spec = catalog->lookup("Ar");
    std::string key = spec->name() + "-" + spec->name();
    pair_potential_ptr_t ptr(new LJ(forceField, box, bc));
    auto pair = std::make_pair(key, ptr);
    map.emplace(pair);

    pair_potential &pairPotential = *map.at(key);
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
    auto forceField = factory::obtainFrom(stream, catalog);
    stream.close();

    test1(catalog, forceField);

    return EXIT_SUCCESS;
}