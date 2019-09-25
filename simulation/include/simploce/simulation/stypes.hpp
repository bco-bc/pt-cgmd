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
 * File:   stypes.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:39 PM
 */

#ifndef STYPES_HPP
#define STYPES_HPP

#include "simploce/particle/ptypes.hpp"
#include "simploce/util/cube.hpp"
#include "simploce/util/param.hpp"
#include <tuple>

namespace simploce {
    
    // Forward declarations.
    class BoundaryCondition;
    class CoarseGrainedForceField;
    class CoarseGrainedDisplacer;
    class AtomisticDisplacer;
    class AtomisticForceField;
    class ModelFactory;
    class ProtonTransferPairListGenerator;
    
    // P is particle type, e.g. Atom or Bead.
    template <typename P>
    class ParticlePairListGenerator;
    
    // P is particle type, e.g. Atom or Bead.
    template <typename P>
    class Interactor;
    
    // P is particle type, e.g. Atom or Bead.
    template <typename P>
    class SimulationModel;    
    
    using temperature_t = simploce::value_t<real_t, 100>;
    
    using density_t = simploce::value_t<real_t, 101>;
    
    using length_t = simploce::value_t<real_t, 102>;
    
    using volume_t = simploce::value_t<real_t, 103>;
    
    using number_density_t = simploce::value_t<real_t, 104>;
    
    using molarity_t = simploce::value_t<real_t, 105>;
    
    using box_t = Cube<real_t>;
    
    using box_ptr_t = std::shared_ptr<box_t>;
    
    using bc_t = BoundaryCondition;
    
    using bc_ptr_t = std::shared_ptr<bc_t>;
    
    using cg_displacer_ptr_t = std::shared_ptr<CoarseGrainedDisplacer>;
    
    using at_displacer_ptr_t = std::shared_ptr<AtomisticDisplacer>;
    
    using at_ff_ptr_t = std::shared_ptr<AtomisticForceField>;
    
    using cg_ff_ptr_t = std::shared_ptr<CoarseGrainedForceField>;
    
    using at_interactor_t = Interactor<Atom>;
        
    using at_interactor_ptr_t = std::shared_ptr<at_interactor_t>;

    using cg_interactor_t = Interactor<Bead>;
        
    using cg_interactor_ptr_t = std::shared_ptr<cg_interactor_t>;
    
    /**
     * Bead pair lists generator type.
     */
    using cg_ppair_list_gen_t = ParticlePairListGenerator<Bead>;
    
    /**
     * Bead pair lists generator pointer type.
     */
    using cg_ppair_list_gen_ptr_t = std::shared_ptr<cg_ppair_list_gen_t>;
    
    /**
     * Atom pair list generator type.
     */
    using at_ppair_list_gen_t = ParticlePairListGenerator<Atom>;
    
    /**
     * Atom pair lists generator pointer type.
     */
    using at_ppair_list_gen_ptr_t = std::shared_ptr<at_ppair_list_gen_t>;
    
    /**
     * Coarse grained simulation model type.
     */
    using cg_sim_model_t = SimulationModel<Bead>;
    
    /**
     * Coarse grained simulation model.
     */
    using cg_sim_model_ptr_t = std::shared_ptr<SimulationModel<Bead>>;
    
    using sim_model_factory_ptr_t = std::shared_ptr<ModelFactory>;
            
    /**
     * Type for holding simulation parameters, such as values for the time step 
     * and the requested number of steps.
     */
    using sim_param_t = param::param_t;
    
    /**
     * Defines a range of something in a given indexed collection to forms pairs between 
     * elements of the collection.
     * tuple_t t;
     * t.get<0> is the begin index of an indexed collection.
     * t.get<1> is the end index.
     * t.get<2> is the number of pairs.
     */
    using range_t = std::tuple<std::size_t, std::size_t, std::size_t>;
    
    /**
     * Protonatable beads pair list generator pointer type.
     */
    using pt_pair_list_gen_ptr_t = std::shared_ptr<ProtonTransferPairListGenerator>;
    
}

#endif /* STYPES_HPP */

