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
    class Particle;
    class Atomistic;
    class CoarseGrained;
    class ParticleSystem;
    class ParticleSpec;
    class ParticleSpecCatalog;
    class ParticleSystemFactory;
    class ParticleGroup;
    class Bond;


    // Type definitions.

    /**
     * Particle pointer type.
     */
     using p_ptr_t = std::shared_ptr<Particle>;

    /**
     * Particle group pointer type.
     */
    using pg_ptr_t = std::shared_ptr<ParticleGroup>;

    /**
     * Atom pointer type.
     */
    using atom_ptr_t = std::shared_ptr<Atom>;

    /**
     * Bead pointer type.
     */
    using bead_ptr_t = std::shared_ptr<Bead>;

    /**
     * Particle system pointer type.
     */
    using p_system_ptr_t = std::shared_ptr<ParticleSystem>;

    /**
     * Atomistic particle system pointer type.
     */
    using at_sys_ptr_t = std::shared_ptr<Atomistic>;

    /**
     * Coarse grained particle system pointer type.
     */
    using cg_sys_ptr_t = std::shared_ptr<CoarseGrained>;

    /**
     * Particle specification pointer type.
     */
    using spec_ptr_t = std::shared_ptr<ParticleSpec>;

    /**
     * Particle specification catalog pointer type.
     */
    using spec_catalog_ptr_t = std::shared_ptr<ParticleSpecCatalog>;

    /**
     * Particle box type.
     */
    using box_t = Cube<real_t>;

    /**
     * Particle box pointer type.
     */
    using box_ptr_t = std::shared_ptr<box_t>;

    /**
     * Pair of particle identifiers type.
     */
    using id_pair_t = std::pair<id_t, id_t>;

    /**
     * Particle system factory pointer type.
     */
    using p_system_fact_ptr_t = std::shared_ptr<ParticleSystemFactory>;





    // OLD

    class ProtonationSiteCatalog;

    template <typename P>
    class ProtonationSite;

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

}

#endif /* P_TYPES_HPP */

