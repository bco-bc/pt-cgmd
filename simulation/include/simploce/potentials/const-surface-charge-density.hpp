/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/29/21.
 */

#ifndef SIMULATION_CONST_SURFACE_CHARGE_DENSITY_HPP
#define SIMULATION_CONST_SURFACE_CHARGE_DENSITY_HPP

#include "external-potential.hpp"
#include "../simulation/s-types.hpp"

namespace simploce {

    /**
     * Energy and force due to a surface carrying a constant uniform surface charge density.
     */
    class ConstantSurfaceChargeDensity : public external_potential {
    public:

        /**
         * Surface location.
         * xy : Surface is in the xy plane at z = 0.
         * yz : Surface is in the yz plane at x = 0;
         * zx : Surface is in the zx plane at y = 0;
         */
        enum PLANE {xy=1, yz, zx};

        /**
         * Conversion from std::string to PLANE.
         * @param value One of "xy", "yz", and "zx".
         * @return PLANE.
         */
        static PLANE valueOf(std::string value);

        ConstantSurfaceChargeDensity(srf_charge_density_t sigma,
                                     real_t eps_r,
                                     PLANE plane,
                                     bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

    private:

        friend class ElectricPotentialDifference;

        static std::pair<energy_t, force_t> forceAndEnergy(srf_charge_density_t sigma,
                                                           PLANE plane,
                                                           real_t eps_r,
                                                           const bc_ptr_t& bc,
                                                           const position_t& r,
                                                           const charge_t& q);

        srf_charge_density_t sigma_;
        real_t eps_r_;
        PLANE plane_;
        bc_ptr_t bc_;
    };
}

#endif //SIMULATION_CONST_SURFACE_CHARGE_DENSITY_HPP
