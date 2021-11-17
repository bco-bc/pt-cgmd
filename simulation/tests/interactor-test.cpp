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

using namespace simploce;

void test1(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    std::cout << "Argon:" << std::endl;
    auto factory = factory::protonatableParticleModelFactory(catalog);
    box_ptr_t box = factory::box(length_t{3.47786});
    auto atomistic = factory->argon(box);

    auto simulationParameters = factory::simulationParameters();
    param::write(std::cout, *simulationParameters);
    auto bc = factory::pbc(atomistic->box());
    auto pairListGenerator = factory::pairListsGeneratorForAtoms(box, bc);

    Interactor<Atom> interactor(simulationParameters, forceField, pairListGenerator, box, bc);
    auto result = interactor.interact(atomistic);
    std::cout << "Non-bonded potential energy: " << result.first << std::endl;
    std::cout << "Bonded potential energy: " << result.second << std::endl;
    std::cout << std::endl;
}

void test2(const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    std::cout << "Diatomic:" << std::endl;
    auto factory = factory::protonatableParticleModelFactory(catalog);
    auto spec = catalog->O();
    auto diatomic = factory->diatomic(0.11, spec);

    auto simulationParameters = factory::simulationParameters();
    param::write(std::cout, *simulationParameters);
    auto box = diatomic->box();
    auto bc = factory::pbc(box);
    auto pairListGenerator = factory::pairListsGeneratorForAtoms(box, bc);

    Interactor<Atom> interactor(simulationParameters, forceField, pairListGenerator, box, bc);
    auto result = interactor.interact(diatomic);
    std::cout << "Non-bonded potential energy: " << result.second << std::endl;
    std::cout << "Bonded potential energy: " << result.first << std::endl;
    std::cout << std::endl;
}

void test3 (const spec_catalog_ptr_t &catalog, const ff_ptr_t &forceField) {
    std::cout << "Polarizable water:" << std::endl;
    auto factory = factory::protonatableParticleModelFactory(catalog);
    auto box = factory::box(7.27);
    auto bc = factory::pbc(box);
    auto polarizableWater = factory->polarizableWater(box);
    auto pairListGenerator = factory::pairListsGeneratorForBeads(box, bc);
    auto simulationParameters = factory::simulationParameters();

    Interactor<Bead> interactor(simulationParameters, forceField, pairListGenerator, box, bc);
    for (int k = 0; k != 100; ++k) {
        std::cout << "Step #" << k << std::endl;
        auto result = interactor.interact(polarizableWater);
        std::cout << "Non-bonded potential energy: " << result.first << std::endl;
        std::cout << "Bonded potential energy: " << result.second << std::endl;
    }
    std::cout << std::endl;
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

    //test1(catalog, forceField);
    //test2(catalog, forceField);
    test3(catalog, forceField);

    return EXIT_SUCCESS;
}
