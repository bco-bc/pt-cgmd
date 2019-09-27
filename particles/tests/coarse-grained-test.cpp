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
#include "simploce/particle/discrete-protonatable-bead.hpp"
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
        std::clog << "Updating..." << std::endl;
        return;
    }
};


void test1()
{
    std::clog << "coarse-grained-test test 1" << std::endl;
    std::clog << "coarse-grained-test test 1" << std::endl;
    
    CoarseGrained cg;
    
    position_t r;
    
    spec_ptr_t spec1_dp = 
        ParticleSpec::create("pspec1", 1.0, 1.0, 1.0);
    spec_ptr_t spec1_p = ParticleSpec::create("dpspec1", 1.0, 1.0, 1.0);
    
    bead_ptr_t bead1 = cg.addBead(1, "test1", r, spec1_dp);
    bead_ptr_t bead2 = cg.addBead(2, "test2", r, spec1_p);
    cg.doWithAll<void>([](const std::vector<bead_ptr_t> &particles) {
        std::clog << "Size: " << particles.size() << std::endl;
        for (auto p : particles) {
            std::clog << " " << p->index() << std::endl;
        }
    });
    
    if ( cg.contains(bead1->index()) ) {
        std::clog << "Particle found: " << *(cg.find(bead1->index())) << std::endl;
    } else {
        std::clog << bead1->index() << ": Not in particle model." << std::endl;
    }
    
    Bond<Bead> bond1 = Bond<Bead>::makeBond(bead1, bead2);
    
    spec_ptr_t pspec = 
        ParticleSpec::create("pspec", 2.0, 1.0, 1.0, 7.0, false);
    dprot_bead_ptr_t pbead = cg.addDiscreteProtonatableBead(3, "ptest", r, 0, pspec);
    
    Bond<Bead> bond2 = Bond<Bead>::makeBond(bead2, pbead);
    
    std::clog << "Bond1 contains bead2? " << bond1.contains(bead2) << std::endl;
    std::clog << "Bond2 contains bead1? " << bond2.contains(bead1) << std::endl;
    
    std::clog << "Size: " << cg.numberOfParticles() << std::endl;
    std::clog << "Charge: " << cg.charge() << std::endl;
    std::clog << "Protonation state: " << cg.protonationState() << std::endl;    
    
    // List specification names.
    cg.doWithAll<void>([](const std::vector<bead_ptr_t> &particles) {
        std::clog << "Particle specification names:" << std::endl;
        for (auto p : particles) {            
            std::clog << p->index() << " " << p->spec()->name() << std::endl;
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
        std::clog << "Particle positions: " << std::endl;
        for (auto p: particles) {
            std::clog << p->index() << " " << p->position() << std::endl;
        }
    });
    
    std::clog << std::endl;
    std::clog << "Reading particle model from input file..." << std::endl;
    std::clog << "Particle specifications:" << std::endl;
    std::ifstream stream;
    file::open_input(stream, 
                    "/home/ajuffer/simploce/pt-cgmd/particles/resources/particles-specs.dat");
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::create(stream);
    stream.close();
    std::clog << *catalog << std::endl;
    
    std::clog << "Particle model:" << std::endl;
    file::open_input(stream,
                    "/home/ajuffer/simploce/pt-cgmd/particles/resources/coarse-grained-system.dat");
    cg_ptr_t cg2 = CoarseGrained::readFrom(stream, catalog);
    stream.close();
    std::clog << "Number of beads: " << cg2->numberOfParticles() << std::endl; 
    std::clog << "State:" << std::endl;
    cg2->writeState(std::clog);
    std::clog << std::endl;
    std::clog << "Particle model:" << std::endl;
    cg2->write(std::clog);
    std::clog << std::endl;

    std::clog << std::endl;
    std::clog << "Update all.." << std::endl;
    cg2->doWithAll<void, Update<Bead>>(Update<Bead>{});
    std::clog << std::endl;
    
    std::clog << std::endl;
    std::clog << "Particle model 'pm':" << std::endl;
    ParticleModel<Bead, ParticleGroup<Bead>>& pm = *cg2;
    pm.doWithAll<void>([] (const std::vector<bead_ptr_t> &particles) {
        std::clog << "Number of particles: " << particles.size() << std::endl;
    });
    
    pm.doWithAllFreeGroups<void>([] (const std::vector<bead_ptr_t>& all,
                                     const std::vector<bead_ptr_t>& free,
                                     const std::vector<bead_group_ptr_t>& groups) {
        std::clog << "Number of particle groups: " << groups.size() << std::endl;
        for (auto g : groups) {
            std::clog << "Group position: " << g->position() << std::endl;
        }
    });
    
    std::clog << "State:" << std::endl;
    cg2->writeState(std::clog);
    std::clog << std::endl;
}


int main(int argc, char** argv) {
    std::clog << "%SUITE_STARTING% coarse-grained-test" << std::endl;
    std::clog << "%SUITE_STARTED%" << std::endl;

    std::clog << "%TEST_STARTED% test1 (coarse-grained-test)" << std::endl;
    test1();
    std::clog << "%TEST_FINISHED% time=0 test1 (coarse-grained-test)" << std::endl;

    std::clog << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

