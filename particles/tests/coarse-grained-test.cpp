/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   coarse-grained-test.cpp
 * Author: ajuffer
 *
 * Created on September 18, 2019, 3:44 PM
 */

#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/protonatable-bead.hpp"
#include "simploce/particle/protonation-site.hpp"
#include "simploce/particle/bond.hpp"
#include "simploce/particle/ptypes.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>

/*
 * Simple C++ Test Suite
 */

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
    std::cout << "coarse-grained-test test 1" << std::endl;
    std::cout << "coarse-grained-test test 1" << std::endl;
    
    CoarseGrained cg;
    
    position_t r;
    
    bead_spec_ptr_t spec1_dp = 
            ParticleSpec::createForBead("pspec1", 1.0, 1.0, 1.0);
    bead_spec_ptr_t spec1_p = ParticleSpec::createForBead("dpspec1", 1.0, 1.0, 1.0);
    
    bead_ptr_t bead1 = cg.addBead("test1", r, spec1_dp);
    bead_ptr_t bead2 = cg.addBead("test2", r, spec1_p);
    cg.doWithAll<void>([](const std::vector<bead_ptr_t> &particles) {
        std::cout << "Size: " << particles.size() << std::endl;
        for (auto p : particles) {
            std::cout << " " << p->id() << std::endl;
        }
    });
    
    if ( cg.contains(bead1->id()) ) {
        std::cout << "Particle found: " << *(cg.find(bead1->id())) << std::endl;
    } else {
        std::cout << bead1->id() << ": Not in particle model." << std::endl;
    }
    
    Bond<Bead> bond1 = Bond<Bead>::makeBond(bead1, bead2);
    
    bead_spec_ptr_t pspec = 
            ParticleSpec::createForProtonatableBead("pspec", 2.0, 1.0, 1.0, 7.0);
    prot_bead_ptr_t pbead = cg.addProtonatableBead("ptest", r, 0, pspec);
    
    Bond<Bead> bond2 = Bond<Bead>::makeBond(bead2, pbead);
    
    std::cout << "Bond1 contains bead2? " << bond1.contains(bead2) << std::endl;
    std::cout << "Bond2 contains bead1? " << bond2.contains(bead1) << std::endl;
    
    std::cout << "Size: " << cg.numberOfParticles() << std::endl;
    std::cout << "Charge: " << cg.charge() << std::endl;
    std::cout << "Protonation state: " << cg.protonationState() << std::endl;    
    
    // List specification names.
    cg.doWithAll<void>([](const std::vector<bead_ptr_t> &particles) {
        std::cout << "Particle specification names:" << std::endl;
        for (auto p : particles) {            
            std::cout << p->id() << " " << p->spec()->name() << std::endl;
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
            std::cout << p->id() << " " << p->position() << std::endl;
        }
    });
    
    std::ifstream stream;
    file::open_input(stream, 
                    "/home/ajuffer/simploce/particles/resources/particles-specs.dat");
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::create(stream);
    stream.close();
    std::clog << *catalog << std::endl;
    
    file::open_input(stream,
                    "/home/ajuffer/simploce/particles/resources/coarse-grained-system.dat");
    cg_ptr_t cg2 = CoarseGrained::createFrom(stream, catalog);
    stream.close();
    std::cout << "Number of beads: " << cg2->numberOfParticles() << std::endl;    
    
    cg2->doWithAll<void, Update<Bead>>(Update<Bead>{});
    
    ParticleModel<Bead, ParticleGroup<Bead>>& pm = *cg2;
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
}

void test2() {
    std::cout << "coarse-grained-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (coarse-grained-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% coarse-grained-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (coarse-grained-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (coarse-grained-test)" << std::endl;

    /*
    std::cout << "%TEST_STARTED% test2 (coarse-grained-test)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (coarse-grained-test)" << std::endl;
     */
    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

