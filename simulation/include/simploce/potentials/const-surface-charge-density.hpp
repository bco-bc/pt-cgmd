/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/29/21.
 */

#ifndef SIMULATION_CONST_SURFACE_CHARGE_DENSITY_HPP
#define SIMULATION_CONST_SURFACE_CHARGE_DENSITY_HPP

#include "external-potential.hpp"
#include "simploce/types/s-types.hpp"
#include "simploce/util/flat-surface.hpp"

namespace simploce {

    /**
     * Energy and force due to a surface carrying a constant uniform surface charge density.
     */
    class ConstantSurfaceChargeDensity : public external_potential {
    public:

        /**
         * Constructor.
         * @param sigma Surface charge density (e/nm^2)
         * @param eps_r Relative permittivity.
         * @param bc Boundary condition.
         * @param flatSurface Flat surface specification. Default is a surface parallel to xy-plane at z = 0.
         */
        ConstantSurfaceChargeDensity(srf_charge_density_t sigma,
                                     FlatSurface flatSurface,
                                     real_t eps_r,
                                     bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

    private:

        friend class ElectricPotentialDifference;

        static std::pair<energy_t, force_t> forceAndEnergy(srf_charge_density_t sigma,
                                                           const FlatSurface& flatSurface,
                                                           real_t eps_r,
                                                           const bc_ptr_t& bc,
                                                           const position_t& r,
                                                           const charge_t& q);

        srf_charge_density_t sigma_;
        real_t eps_r_;
        bc_ptr_t bc_;
        FlatSurface flatSurface_;
    };
}

#endif //SIMULATION_CONST_SURFACE_CHARGE_DENSITY_HPP
