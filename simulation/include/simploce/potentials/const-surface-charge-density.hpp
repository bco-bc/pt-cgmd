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
         * @param eps_r Relative permittivity, screens the electric field.
         * @param bc Boundary condition.
         * @param flatSurface Flat surface specification.
         * @param mesoscopic If true, this potential is for mesoscopic simulations (e.g., DPD).
         */
        ConstantSurfaceChargeDensity(srf_charge_density_t sigma,
                                     FlatSurface flatSurface,
                                     real_t eps_r,
                                     bc_ptr_t bc,
                                     bool mesoscopic = false);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

        /**
         * Returns energy and force on particle.
         * @param ro Position, possible outside simulation box.
         * @param q Charge value.
         * @param df Damping factor for electrostatic interactions when dealing with charge densities. The
         * default value implies the charge is taken as an ideal point charge.
         * @return
         */
        std::pair<energy_t, force_t> forceAndEnergy(const position_t& ro, const charge_t& q, real_t df = 0);

    private:

        srf_charge_density_t sigma_;
        FlatSurface flatSurface_;
        real_t eps_r_;
        bc_ptr_t bc_;
        bool mesoscopic_;
    };
}

#endif //SIMULATION_CONST_SURFACE_CHARGE_DENSITY_HPP
