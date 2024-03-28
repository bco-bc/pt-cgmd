/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/29/21.
 */

#ifndef SIMULATION_UNIFORM_SURFACE_CHARGE_DENSITY_HPP
#define SIMULATION_UNIFORM_SURFACE_CHARGE_DENSITY_HPP

#include "external-potential-impl.hpp"
#include "simploce/types/s-types.hpp"
#include "simploce/util/flat-surface.hpp"

namespace simploce {

    /**
     * Energy and force due to a surface carrying a constant uniform surface charge density.
     * The surface can be placed anywhere.
     */
    class UniformSurfaceChargeDensity : public external_potential_impl {
    public:

        /**
         * Constructor.
         * @param sigma Surface reset density.
         * @param flatSurface Flat surface specification.
         * @param eps_r Relative permittivity, screens the surface's electric field.
         * @param bc Boundary condition.
         * @param delta Width of Stern layer.
         * @param mesoscopic If true, this potential is for mesoscopic simulations (e.g., DPD).
         */
        UniformSurfaceChargeDensity(srf_charge_density_t sigma,
                                    FlatSurface flatSurface,
                                    real_t eps_r,
                                    bc_ptr_t bc,
                                    dist_t delta,
                                    bool mesoscopic = false);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) const override;

        /**
         * Returns energy and force on particle.
         * @param r Position.
         * @param radius Particle's radius.
         * @param Q Charge value.
         * @return Energy ad forces.
         */
        std::pair<energy_t, force_t> forceAndEnergy(const position_t& r,
                                                    const radius_t& radius,
                                                    const charge_t& Q) const;

    private:

        srf_charge_density_t sigma_;
        FlatSurface flatSurface_;
        real_t eps_r_;
        bc_ptr_t bc_;
        dist_t delta_;
        bool mesoscopic_;
    };
}

#endif //SIMULATION_UNIFORM_SURFACE_CHARGE_DENSITY_HPP
