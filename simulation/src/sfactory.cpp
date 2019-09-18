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
 * File:   factory.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 3:41 PM
 */

#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/cg-pol-water.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/simulation/stypes.hpp"
#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/model-factory.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/pfactory.hpp"
#include <memory>

namespace simploce {
    namespace factory {
        
        using at_leap_frog_t = LeapFrog<Atomistic>;
        using cg_leap_frog_t = LeapFrog<CoarseGrained>;
        using at_vv_t = VelocityVerlet<Atomistic>;
        using cg_vv_t = VelocityVerlet<CoarseGrained>;
        using at_lvv_t = LangevinVelocityVerlet<Atomistic>;
        using cg_lvv_t = LangevinVelocityVerlet<CoarseGrained>;
        
        // Coarse grained polarizable force field.
        static cg_ff_ptr_t cgPolWaterFF_{};
        
        // Pair lists generators.
        static cg_ppair_list_gen_ptr_t  cgPairlistGen_{};
        static at_ppair_list_gen_ptr_t atpairListGen_{};
        
        // Leap Frog
        static at_displacer_ptr_t atLeapFrog_{};
        static cg_displacer_ptr_t cgLeapFrog_{};
        
        // Velocity Verlet
        static at_displacer_ptr_t atVV_{};
        static cg_displacer_ptr_t cgVV_{};
        
        // Velocity Verlet
        static at_displacer_ptr_t atLVV_{};
        static cg_displacer_ptr_t cgLVV_{};
        
        static sim_model_factory_ptr_t modelFactory_{};
        
        static box_ptr_t box_{};
        static bc_ptr_t pbc_{};
        
        // Interactor for coarse grained polarizable water.
        static std::shared_ptr<Interactor<Bead>> cgPolWaterInteractor_{};
                
        static cg_ff_ptr_t 
        coarseGrainedPolarizableWaterForceField(const spec_catalog_ptr_t& catalog,
                                                const bc_ptr_t& bc)
        {
            if ( !cgPolWaterFF_ ) {
                cgPolWaterFF_ = 
                    std::make_shared<CoarseGrainedPolarizableWater>(catalog, bc);
            }
            return cgPolWaterFF_;
        }
        
        sim_model_factory_ptr_t modelFactory(const spec_catalog_ptr_t& catalog)
        {
            if ( !modelFactory_ ) {
                modelFactory_ = std::make_shared<ModelFactory>(catalog);
            }
            return modelFactory_;
        }
        
        cg_ppair_list_gen_ptr_t coarseGrainedPairListGenerator(const box_ptr_t& box,
                                                               const bc_ptr_t& bc)
        {
            if ( !cgPairlistGen_ ) {
                cgPairlistGen_ = std::make_shared<ParticlePairListGenerator<Bead>>(box, bc);
            }
            return cgPairlistGen_;
        }
        
        at_ppair_list_gen_ptr_t atomisticPairListGenerator(const box_ptr_t& box,
                                                           const bc_ptr_t& bc)
        {
            if ( !atpairListGen_ ) {
                atpairListGen_ = std::make_shared<ParticlePairListGenerator<Atom>>(box, bc);
            }
            return atpairListGen_;
        }
        
        cg_interactor_ptr_t 
        interactorCoarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog,
                                                const box_ptr_t& box, 
                                                const bc_ptr_t& bc)        
        {
            cg_ppair_list_gen_ptr_t generator = 
                    factory::coarseGrainedPairListGenerator(box, bc);
            if ( !cgPolWaterInteractor_ ) {
                cg_ff_ptr_t forcefield = coarseGrainedPolarizableWaterForceField(catalog, bc);
                cgPolWaterInteractor_ = std::make_shared<Interactor<Bead>>(forcefield, generator);
            }
            return cgPolWaterInteractor_;
        }
        
        at_displacer_ptr_t leapFrog(at_interactor_ptr_t& interactor)
        {
            if ( !atLeapFrog_ ) {
                atLeapFrog_ = std::make_shared<at_leap_frog_t>(interactor);
            }
            return atLeapFrog_;
        }
        
        cg_displacer_ptr_t leapFrog(cg_interactor_ptr_t& interactor)
        {
            if ( !cgLeapFrog_ ) {
                cgLeapFrog_ = std::make_shared<cg_leap_frog_t>(interactor);
            }
            return cgLeapFrog_;
            
        }
        
        at_displacer_ptr_t velocityVerlet(at_interactor_ptr_t& interactor)
        {
            if ( !atVV_ ) {
                atVV_ = std::make_shared<at_vv_t>(interactor);
            }
            return atVV_;            
        }
        
        cg_displacer_ptr_t velocityVerlet(cg_interactor_ptr_t& interactor)
        {
            if ( !cgVV_ ) {
                cgVV_ = std::make_shared<cg_vv_t>(interactor);
            }
            return cgVV_;                        
        }
        
        at_displacer_ptr_t langevinVelocityVerlet(at_interactor_ptr_t& interactor)
        {
            if ( !atLVV_ ) {
                atVV_ = std::make_shared<at_lvv_t>(interactor);
            }
            return atLVV_;                        
        }
        
        cg_displacer_ptr_t langevinVelocityVerlet(cg_interactor_ptr_t& interactor)
        {
            if ( !cgLVV_ ) {
                cgLVV_ = std::make_shared<cg_lvv_t>(interactor);
            }
            return cgLVV_;                                    
        }
        
        box_ptr_t cube(const length_t& side)
        {
            if ( !box_ ) {
                box_ = std::make_shared<Cube<real_t>>(side());
            }
            return box_;
        }
        
        bc_ptr_t pbc(const box_ptr_t& box)
        {
            if ( !pbc_ ) {
                pbc_ = std::make_shared<PeriodicBoundaryCondition>(box);
            }
            return pbc_;
        }
        
        at_ptr_t atomistic()
        {
            return std::make_shared<Atomistic>();
        }
        
        cg_ptr_t coarseGrained()
        {
            return std::make_shared<CoarseGrained>();
        }
        
    }
}
