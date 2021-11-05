/*
 * File:   coarse-grained-argon.cpp
 * Author: André H. Juffer, Biocenter Oulu
 *
 * Created on September 18, 2019, 3:44 PM
 */

#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/bond.hpp"
#include "simploce/particle/p-types.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

template <typename P>
struct Update {
    void operator () (const std::vector<std::shared_ptr<P>>& particles) const
    {
        std::cout << "Updating..." << std::endl;
        return;
    }
};

void test1()
{
    CoarseGrained cg;
    
    position_t r;

    spec_ptr_t spec1_dp =
        ParticleSpec::create("pspec1", 1.0, 1.0, 1.0, "dp");
    spec_ptr_t spec1_p = ParticleSpec::create("dpspec1", 1.0, 1.0, 1.0, "p");
    
    bead_ptr_t bead1 = cg.addBead("argon", spec1_dp);
    bead1->position(r);
    bead_ptr_t bead2 = cg.addBead("test2", spec1_p);
    bead2->position(r);
    cg.doWithAll<void>([](const std::vector<bead_ptr_t> &particles) {
        std::cout << "Size: " << particles.size() << std::endl;
        for (auto p : particles) {
            std::cout << " " << p->index() << std::endl;
        }
    });
    
    if ( cg.contains(bead1->id()) ) {
        std::cout << "Particle found: " << *(cg.find(bead1->id())) << std::endl;
    } else {
        std::cout << bead1->index() << ": Not in particle model." << std::endl;
    }
    
    Bond<Bead> bond1 = Bond<Bead>::makeBond(bead1, bead2);
    
    spec_ptr_t pspec = 
        ParticleSpec::create("pspec", 2.0, 1.0, 1.0, 7.0, "argon");
    auto pbead = cg.addBead("ptest", pspec);
    pbead->position(r);
    
    Bond<Bead> bond2 = Bond<Bead>::makeBond(bead2, pbead);
    
    std::cout << "Bond1 contains bead2? " << bond1.contains(bead2) << std::endl;
    std::cout << "Bond2 contains bead1? " << bond2.contains(bead1) << std::endl;
    
    std::cout << "Size: " << cg.numberOfParticles() << std::endl;
    std::cout << "Charge: " << cg.charge() << std::endl;
    
    // List specification names.
    cg.doWithAll<void>([](const std::vector<bead_ptr_t> &particles) {
        std::cout << "Particle specification names:" << std::endl;
        for (auto p : particles) {            
            std::cout << p->index() << " " << p->spec()->name() << std::endl;
        }

    });
    
    cg.doWithAll<void>([](const std::vector<bead_ptr_t> &particles) {
        position_t r{0.0, 1.0, 2.0};
        for (auto p : particles) {
            r *= 2.0;
            p->position(r);
        }
    });
    
    cg.doWithAll<void>([] (const std::vector<bead_ptr_t> &particles) {
        std::cout << "Particle positions: " << std::endl;
        for (auto p: particles) {
            std::cout << p->index() << " " << p->position() << std::endl;
        }
    });
    
    std::cout << std::endl;
    std::cout << "Reading particle model from input file..." << std::endl;
    std::cout << "Particle specifications:" << std::endl;
    std::ifstream stream;
    util::open_input_file(stream,
                    "/localdisk/resources/particles-specs.dat");
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::obtainFrom(stream);
    stream.close();
    std::cout << *catalog << std::endl;
    
    std::cout << "Particle model:" << std::endl;
    util::open_input_file(stream,
                    "/localdisk/resources/coarse-grained-system.dat");
    cg_mod_ptr_t cg2 = CoarseGrained::obtainFrom(stream, catalog);
    stream.close();
    std::cout << "Number of beads: " << cg2->numberOfParticles() << std::endl; 
    std::cout << "State:" << std::endl;
    cg2->writeState(std::cout);
    std::cout << std::endl;
    std::cout << "Particle model:" << std::endl;
    cg2->write(std::cout);
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "Update all.." << std::endl;
    cg2->doWithAll<void, Update<Bead>>(Update<Bead>{});
    std::cout << std::endl;
    
    std::cout << std::endl;
    std::cout << "Particle model 'pm':" << std::endl;
    ParticleSystem<Bead, ParticleGroup<Bead>>& pm = *cg2;
    pm.doWithAll<void>([] (const std::vector<bead_ptr_t> &particles) {
        std::cout << "Number of particles: " << particles.size() << std::endl;
    });
    
    pm.doWithAllFreeGroups<void>([] (const std::vector<bead_ptr_t>& all,
                                     const std::vector<bead_ptr_t>& free,
                                     const std::vector<bead_group_ptr_t>& groups) {
        std::cout << "Number of particle groups: " << groups.size() << std::endl;
        for (auto g : groups) {
            std::cout << "Group position: " << g->position() << std::endl;
        }
    });
    
    std::cout << "State:" << std::endl;
    cg2->writeState(std::cout);
    std::cout << std::endl << std::endl;
    std::cout << "Empty?: " << cg2->empty() << std::endl;

    cg2->resetForces();
}


int main() {
    try {
        test1();
    } catch (std::exception &exception) {
        std::cerr << "Test halted due to exception. " << exception.what() << std::endl;
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

