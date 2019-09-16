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
 * Defines all types in use by Particles.
 * File:   types.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 11:54 AM
 */

#ifndef TYPES_HPP
#define TYPES_HPP

#include "simploce/util/cvector_t.hpp"
#include "simploce/util/value_t.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    // Forward declarations.
    class Atom;
    class Bead;
    class ProtonatableBead;
    class ParticleSpec;    
    class ParticleSpecCatalog;
    class ProtonationSiteCatalog;
    class Atomistic;
    class CoarseGrained;
    
    template <typename T>
    class ParticleGroup;
    
    template <typename T, typename S>
    class ProtonationSite;
    
    /**
     * Real number.
     */
    using real_t = double;
    
    /**
     * Charge value.
     */
    using charge_t = simploce::value_t<real_t, 1>;
    
    /**
     * Mass value.
     */
    using mass_t = simploce::value_t<real_t, 2>;
    
    /**
     * pKa value.
     */
    using pKa_t = simploce::value_t<real_t, 3>;
    
    /**
     * Radius.
     */
    using radius_t = simploce::value_t<real_t, 4>;
    
    /**
     * Time.
     */
    using stime_t = simploce::value_t<real_t, 5>;
    
    /**
     * Energy.
     */
    using energy_t = simploce::value_t<real_t, 6>;
    
    /**
     * Position.
     */
    using position_t = cvector_t<real_t, 1>;
    
    /**
     * Linear momentum.
     */
    using momentum_t = cvector_t<real_t, 2>;
    
    /**
     * Force
     */
    using force_t = cvector_t<real_t, 3>;
    
    /**
     * Distance (difference) type
     */
    using dist_vect_t = position_t;        
    
    /**
     * Atomistic model pointer type.
     */
    using at_ptr_t = std::shared_ptr<Atomistic>;
    
    /**
     * Atom pointer type
     */
    using atom_ptr_t = std::shared_ptr<Atom>;
    
    /**
     * Coarse grained model pointer type.
     */
    using cg_ptr_t = std::shared_ptr<CoarseGrained>;
    
    /**
     * Bead pointer type.
     */
    using bead_ptr_t = std::shared_ptr<Bead>;        
    
    /**
     * Protonatable bead pointer type.
     */
    using prot_bead_ptr_t = std::shared_ptr<ProtonatableBead>;
    
    /**
     * Particle specification pointer type
     */
    using particle_spec_ptr_t = std::shared_ptr<ParticleSpec>;     
    
    /**
     * Atom specification type.
     */
    using AtomSpec = ParticleSpec;
    
    /**
     * Atom specification type.
     */
    using atom_spec_ptr_t = particle_spec_ptr_t;
    
    /**
     * Bead specification type.
     */
    using BeadSpec = ParticleSpec;
    
    /**
     * Bead specification pointer type.
     */
    using bead_spec_ptr_t = particle_spec_ptr_t;
    
    /**
     * Protonatable bead specification.
     */
    using ProtBeadSpec = ParticleSpec;
    
    /**
     * Protonatable bead specification pointer type.
     */
    using prot_bead_spec_ptr_t = particle_spec_ptr_t;
    
    /**
     * Particle specification catalog pointer type.
     */
    using spec_catalog_ptr_t = std::shared_ptr<ParticleSpecCatalog>;
    
    /**
     * Atom-based protonation site type.
     */
    using atom_prot_site_t = ProtonationSite<Atom, AtomSpec>;
    
    /**
     * Atom-based protonation site pointer type.
     */    
    using atom_prot_site_ptr_t = std::shared_ptr<atom_prot_site_t>;
    
    /**
     * Bead-based protonation site type.
     */
    using bead_prot_site_t = ProtonationSite<Bead, BeadSpec>;
    
    /**
     * Bead-based protonation site pointer type.
     */    
    using bead_prot_site_ptr_t = std::shared_ptr<bead_prot_site_t>;
    
    /**
     * Protonation sites catalog pointer type.
     */
    using prot_site_catalog_ptr_t = std::shared_ptr<ProtonationSiteCatalog>;
    
    /**
     * Atom group type.
     */
    using atom_group_t = ParticleGroup<Atom>;
            
    /**
     * Atom group pointer type.
     */
    using atom_group_ptr_t = std::shared_ptr<atom_group_t>;
    
    /**
     * Bead group type.
     */
    using bead_group_t = ParticleGroup<Bead>;
            
    /**
     * Atom group pointer type.
     */
    using bead_group_ptr_t = std::shared_ptr<bead_group_t>;
    
    /**
     * Type for a pair of particle identifiers.
     */
    using id_pair_t = std::pair<int, int>;
    
}

#endif /* TYPES_HPP */

