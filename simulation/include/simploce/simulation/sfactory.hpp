/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   factory.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 3:39 PM
 */

#ifndef SFACTORY_HPP
#define SFACTORY_HPP

#include "stypes.hpp"
#include "simploce/particle/pfactory.hpp"

namespace simploce {
    namespace factory {
        
        /**
         * Coarse grained force field for polarizable water.
         * @param catalog Particle specification catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @param protonatable If true, water is protonatable.
         * @return Force field.
         */
        cg_ff_ptr_t
        polarizableWaterForceField(const spec_catalog_ptr_t& catalog,
                                   const bc_ptr_t& bc,
                                   const box_ptr_t& box,
                                   bool protonatable = false);
        
        /**
         * Returns coarse grained force field for formic acid (HCOOH) in water.
         * @param catalog Particle specification catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @param water Coarse grained force field for water.
         * @return Force field.
         */
        cg_ff_ptr_t 
        formicAcidSolutionForceField(const spec_catalog_ptr_t& catalog,
                                     const bc_ptr_t& bc,
                                     const box_ptr_t& box,
                                     const cg_ff_ptr_t& water);
        
        /**
         * Returns coarse grained force field for an NaCl electrolyte solution.
         * @param catalog Particle specification catalog.
         * @param bc Boundary condition.
         * @param box Simulation box.
         * @return Force field.
         */
        cg_ff_ptr_t
        electrolyteForceField(const spec_catalog_ptr_t& catalog,
                              const bc_ptr_t& bc,
                              const box_ptr_t& box);
        
        /**
         * Returns simulation model factory.
         * @param particle_model_fact_ptr_t Particle model factory.
         * @param catalog Particle specifications catalog.
         * @return Model factory.
         */
        sim_model_fact_ptr_t 
        simulationModelFactory(const spec_catalog_ptr_t& catalog);
        
        /**
         * Returns particle pair list generator for coarse grained particle models.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Pair list generator.
         */
        cg_ppair_list_gen_ptr_t
        coarseGrainedPairListGenerator(const box_ptr_t& box,
                                       const bc_ptr_t& bc);
        
        /**
         * Returns particle pair list generator for atomistic particle models.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Pair list generator.
         */
        at_ppair_list_gen_ptr_t 
        atomisticPairListGenerator(const box_ptr_t& box,
                                   const bc_ptr_t& bc);
        
        
        /**
         * Returns coarse grained interactor.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @param protonatable If true, the interactor will assume parameters for 
         * protonatable waters.
         * @return Interactor.
         */
        cg_interactor_ptr_t 
        polarizableWaterInteractor(const spec_catalog_ptr_t& catalog,
                                   const box_ptr_t& box, 
                                   const bc_ptr_t& bc,
                                   bool protonatable = false);
        
        /**
         * Returns coarse grained interactor.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @return Interactor.
         */
        cg_interactor_ptr_t
        formicAcidSolutionInteractor(const spec_catalog_ptr_t& catalog,
                                     const box_ptr_t& box, 
                                     const bc_ptr_t& bc);
        
        /**
         * Returns coarse grained interactor for electrolyte solution.
         * @param catalog Particle specifications catalog.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @return Interactor.
         */
        cg_interactor_ptr_t
        electrolyteInteractor(const spec_catalog_ptr_t& catalog,
                              const box_ptr_t& box, 
                              const bc_ptr_t& bc);
        
        
        /**
         * Leap frog algorithm for atomistic particle models.
         * @param interactor Atomistic interactor.
         * @return Algorithm.
         */
        at_displacer_ptr_t 
        leapFrog(at_interactor_ptr_t& interactor);
        
        /**
         * Leap frog algorithm for coarse grained particle models.
         * @param interactor Coarse grained interactor.
         * @return Algorithm.
         */
        cg_displacer_ptr_t 
        leapFrog(cg_interactor_ptr_t& interactor);
        
        /**
         * Velocity Verlet algorithm for atomistic particle models.
         * @param interactor Atomistic interactor.
         * @return Algorithm.
         */
        at_displacer_ptr_t 
        velocityVerlet(at_interactor_ptr_t& interactor);
        
        /**
         * Velocity Verlet algorithm for coarse grained particle models.
         * @param interactor Coarse grained interactor.
         * @return Algorithm.
         */
        cg_displacer_ptr_t 
        velocityVerlet(cg_interactor_ptr_t& interactor);

        /**
         * Langevin Velocity Verlet algorithm for atomistic particle models.
         * @param interactor Atomistic interactor.
         * @return Algorithm.
         */
        at_displacer_ptr_t 
        langevinVelocityVerlet(at_interactor_ptr_t& interactor);
        
        /**
         * Langevin Velocity Verlet algorithm for coarse grained particle models.
         * @param interactor Coarse grained interactor.
         * @return Algorithm.
         */
        cg_displacer_ptr_t 
        langevinVelocityVerlet(cg_interactor_ptr_t& interactor);
        
        /**
         * Returns Langevin Velocity Verlet algorithm with proton transfer for 
         * coarse grained models.
         * @param interactor Coarse grained interactor.
         * @param generator Protonatable bead pair list generator.
         * @param displacer Proton transfer displacer.
         * @return Displacer.
         */
        cg_displacer_ptr_t
        protonTransferlangevinVelocityVerlet(const cg_interactor_ptr_t& interactor,
                                             const pt_pair_list_gen_ptr_t& generator,
                                             const pt_displacer_ptr_t& displacer);
        
        /**
         * Returns periodic boundary conditions.
         * @param box Simulation box.
         * @return Boundary condition.
         */
        bc_ptr_t 
        pbc(const box_ptr_t& box);
        
        /**
         * Returns generator of pairs of protonatable beads possibly involved in proton
         * transfer.
         * @param rmax Max distance between protonatable beads.
         * @param bc Boundary condition.
         * @return Generator.
         */
        pt_pair_list_gen_ptr_t 
        protonTransferPairListGenerator(const length_t& rmax,
                                        const bc_ptr_t& bc);
        
        /**
         * Returns proton transfer (PT) displacer with constant rate.
         * @param rate Rate.
         * @param gamma Inverse of time constant of decay.
         * @return PT displacer
         */
        pt_displacer_ptr_t constantRate(const rate_t& rate, const real_t& gamma);
    }
}

#endif /* FACTORY_HPP */

