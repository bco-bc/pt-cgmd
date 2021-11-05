/*
 * Various types associated with particles.
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 11:54 AM
 */

#ifndef P_TYPES_HPP
#define P_TYPES_HPP

#include "simploce/types/u-types.hpp"
#include "simploce/util/cube.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    // Forward declarations.

    class Atom;
    class Bead;
    class ParticleSpec;
    class ParticleSpecCatalog;

    template <typename P>
    class ParticleGroup;


    // Type definitions.

    /**
     * Atom pointer type.
     */
    using atom_ptr_t = std::shared_ptr<Atom>;

    /**
     * Bead pointer type.
     */
    using bead_ptr_t = std::shared_ptr<Bead>;

    /**
     * Particle specification pointer type.
     */
    using spec_ptr_t = std::shared_ptr<ParticleSpec>;

    /**
     * Particle specification catalog pointer type.
     */
    using spec_catalog_ptr_t = std::shared_ptr<ParticleSpecCatalog>;



    class ProtonationSiteCatalog;
    class Atomistic;
    class CoarseGrained;
    class PolarizableWater;
    class ParticleModelFactory;


    template <typename P>
    class ProtonationSite;


    /**
     * Particle box type.
     */
    using box_t = Cube<real_t>;
    
    /**
     * Particle box pointer type.
     */
    using box_ptr_t = std::shared_ptr<box_t>;

    /**
     * Atomistic particle model pointer type.
     */
    using at_mod_ptr_t = std::shared_ptr<Atomistic>;

    /**
     * Coarse grained particle model pointer type.
     */
    using cg_mod_ptr_t = std::shared_ptr<CoarseGrained>;
    
    /**
     * Coarse grained polarizable particle model pointer type.
     */
    using cg_pol_water_ptr_t = std::shared_ptr<PolarizableWater>;
    

    /**
     * Atom-based protonation site.
     */
    using atom_prot_site_t = ProtonationSite<Atom>;
    
    /**
     * Atom-based protonation site pointer type.
     */    
    using atom_prot_site_ptr_t = std::shared_ptr<atom_prot_site_t>;
    
    /**
     * Bead-based protonation site type.
     */
    using bead_prot_site_t = ProtonationSite<Bead>;
    
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
     * Bead group pointer type.
     */
    using bead_group_ptr_t = std::shared_ptr<bead_group_t>;
    
    /**
     * Pair of particle identifiers type.
     */
    using id_pair_t = std::pair<id_t, id_t>;
    
    /**
     * Particle model factory pointer type.
     */
    using particle_model_fact_ptr_t = std::shared_ptr<ParticleModelFactory>;
    
}

#endif /* P_TYPES_HPP */

