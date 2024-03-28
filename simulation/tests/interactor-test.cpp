/*
 * Author: André H. Juffer.
 * Created on 10/11/2021, 16:46.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/interactor.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/param.hpp"
#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"

using namespace simploce;

void mc(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    std::cout << "Argon:" << std::endl;
    auto factory = factory::protonatableParticleSystemFactory(catalog);
    box_ptr_t box = factory::box(length_t{3.47786});
    auto atomistic = factory->argon(box);

    auto param = factory::simulationParameters();
    param::write(std::cout, *param);
    auto bc = factory::pbc(atomistic->box());
    auto cutoff = param->get<real_t>("forces.nb.cutoff");
    auto pairListGenerator = factory::pairListGenerator(param, bc);

    auto forces = factory::forces(param, bc, forceField);
    auto interactor = factory::interactor(param, forceField, bc);
    auto result = interactor->interact(atomistic);
    std::cout << "Bonded potential energy: " << std::get<0>(result) << std::endl;
    std::cout << "Non-bonded potential energy: " << std::get<1>(result)<< std::endl;
    std::cout << "External potential energy: " << std::get<2>(result)<< std::endl;
    std::cout << std::endl;
}

void test2(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    std::cout << "Diatomic:" << std::endl;
    auto factory = factory::protonatableParticleSystemFactory(catalog);
    auto spec = catalog->O();
    auto diatomic = factory->diatomic(0.11, spec);

    auto simulationParameters = factory::simulationParameters();
    param::write(std::cout, *simulationParameters);
    auto box = diatomic->box();
    auto bc = factory::pbc(box);

    auto pairListGenerator = factory::pairListGenerator(simulationParameters, bc);

    auto interactor = factory::interactor(simulationParameters, forceField, bc);
    auto result = interactor->interact(diatomic);
    std::cout << "Bonded potential energy: " << std::get<0>(result) << std::endl;
    std::cout << "Non-bonded potential energy: " << std::get<1>(result)<< std::endl;
    std::cout << "External potential energy: " << std::get<2>(result)<< std::endl;
    std::cout << std::endl;
}

void test3 (const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    std::cout << "Polarizable water:" << std::endl;
    auto factory = factory::protonatableParticleSystemFactory(catalog);
    auto box = factory::box(7.27);
    auto bc = factory::pbc(box);
    auto polarizableWater = factory->cgmdPolarizableWater(box);
    auto param = factory::simulationParameters();
    auto cutoff = param->get<real_t>("forces.nb.cutoff");
    auto pairListGenerator = factory::pairListGenerator(param, bc);

    auto interactor = factory::interactor(param, forceField, bc);
    auto result = interactor->interact(polarizableWater);
    std::cout << "Bonded potential energy: " << std::get<0>(result) << std::endl;
    std::cout << "Non-bonded potential energy: " << std::get<1>(result)<< std::endl;
    std::cout << "External potential energy: " << std::get<2>(result)<< std::endl;

    std::cout << std::endl;
}

int main() {
    //util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fnSpecs = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fnSpecs);
    auto catalog = factory::particleSpecCatalog(stream);
    stream.close();

    std::string fnInteractions = "/localdisk/resources/interaction-parameters.dat";
    util::open_input_file(stream, fnInteractions);
    auto forceField = factory::forceField(stream, catalog);
    stream.close();

    //test1(catalog, forceField);
    //test2(catalog, forceField);
    test3(catalog, forceField);

    return EXIT_SUCCESS;
}

