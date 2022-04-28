/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 12/20/21
 */

#ifndef SIMULATION_WALL_HPP
#define SIMULATION_WALL_HPP

#include "external-potential.hpp"
#include "simploce/util/flat-surface.hpp"

namespace simploce {

    /**
     * Impenetrable wall. The associated potential is a simple Lennard-Jones potential. An
     * electrostatic component, for instance, according to a simple Gouy-Chapman theory, may be added.
     */
    class Wall : public external_potential {
    public:

        /**
         * Constructor.
         * @param C12 C12 Lennard-Jones interaction parameter.
         * @param C6 C6 Lennard-Jones interaction parameter.
         * @param flatSurface Surface specification. Default is a surface parallel to the xy plane at z=0.
         * @param sigma Surface (plane) surface charge density. Default is an uncharged surface.
         */
        Wall(real_t C12, real_t C6, bc_ptr_t bc, FlatSurface flatSurface = FlatSurface{}, srf_charge_density_t sigma = 0.0);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

    private:

        real_t C12_;
        real_t C6_;
        bc_ptr_t bc_;
        FlatSurface flatSurface_;
        srf_charge_density_t sigma_;

    };
}

#endif //SIMULATION_WALL_HPP
