/*
 * Author: André H. Juffer.
 * Created on 07/08/2023
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_SOLID_SPHERE_DSF_HPP
#define SIMULATION_SOLID_SPHERE_DSF_HPP

#include "simploce/potentials/pair-potential.hpp"

namespace simploce {

    /**
     * Damped shifted force Coulomb interaction. "Damping" occurs at short distance as well as longer distances.
     * At the cutoff distance and beyond it, the force is exactly zero, while at shorter distances the field
     * remains finite.
     */
    class SolidSphere_DSF : public pair_potential {
    public:

        /**
         * Constructor. All arguments are required.
         * @param cutoff Cutoff distance.
         * @param forceField Force field.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @param radius Radius solid sphere.
         * @param mesoscopic If true, assume mesoscopic units.
         */
        SolidSphere_DSF(dist_t cutoff,
                        ff_ptr_t forceField,
                        box_ptr_t box,
                        bc_ptr_t bc,
                        dist_t radius,
                        bool mesoscopic = false);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& rij,
                                                    real_t Rij,
                                                    real_t Rij2,
                                                    const charge_t& q1,
                                                    const charge_t& q2);

        dist_t cutoff_;
        ff_ptr_t forceField_;
        box_ptr_t box_;
        bc_ptr_t bc_;
        dist_t radius_;
        bool mesoscopic_;
    };
}

#endif //SIMULATION_SOLID_SPHERE_DSF_HPP
