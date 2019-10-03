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
         * @param box SImulation box.
         * @param bc Boundary condition.
         * @return Interactor.
         */
        cg_interactor_ptr_t 
        interactorCoarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog,
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
    }
}

#endif /* FACTORY_HPP */

