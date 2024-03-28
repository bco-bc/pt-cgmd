/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 28 February 2024. Adapted from original version (1996).
 * @see https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(199612)17:16%3C1783::AID-JCC1%3E3.0.CO;2-J
 */

#ifndef SIMULATION_VIRTUAL_PLANE_HPP
#define SIMULATION_VIRTUAL_PLANE_HPP

#include "external-potential-impl.hpp"

namespace simploce {

    /**
     * Holds surface charge density accumulated in the course of a simulation for systems with
     * 2D symmetry in the x- and y-direction, but not in the z-direction, to model asymmetry in
     * charge density outside the central simulation box. For example, an
     * electrolyte solution next to a charged surface parallel the xy-plane at z=0,
     * where a virtual plane would represent the accumulated charge density in a slab of
     * given width placed parallel to the charged surface. The plane is located at specified
     * distance to that charged surface. A particle inside the box will interact with the
     * -average- surface charge density that is -outside- the simulation box.
     */
    class VirtualPlane  {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @param location Distance to the xy-plane at z = 0.
         * @param eps_r Relative permittivity.
         */
        explicit VirtualPlane(box_ptr_t box,
                              bc_ptr_t bc,
                              dist_t location,
                              real_t eps_r);

        /**
         * Returns surface charge density.
         * @return Surface charge density.
         */
        srf_charge_density_t
        surfaceChargeDensity() const;

        /**
         * Returns location.
         * @return Distance to the xy-plane at z=0.
         */
        dist_t
        location() const;

        /**
         * Interaction energy only. Forces not yet provided.
         */
        std::pair<energy_t, force_t>
        operator () (const p_ptr_t& particle) const;

    private:

        friend class VirtualPlanes;

        /**
         * Resets the surface charge density to the given value.
         * @param sigma New value.
         */
        void reset(srf_charge_density_t sigma);

        box_ptr_t box_;
        bc_ptr_t bc_;
        dist_t location_;             // Distance to the xy-plane placed at z=0.
        real_t eps_r_;                // Relative permittivity.
        srf_charge_density_t sigma_;  // Average surface charge density.
    };

    /**
     * Writes virtual surface to output stream, that is location and surface charge density.
     * @param stream Output stream.
     * @param plane Virtual plane.
     * @return Output stream.
     */
    std::ostream&
    operator << (std::ostream& stream, const VirtualPlane& plane);
}


#endif //SIMULATION_VIRTUAL_PLANE_HPP
