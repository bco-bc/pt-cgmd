/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 28 February 2024. Adapted from original version (1996).
 * @see https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(199612)17:16%3C1783::AID-JCC1%3E3.0.CO;2-J
 */


#ifndef SIMULATION_VIRTUAL_PLANES_HPP
#define SIMULATION_VIRTUAL_PLANES_HPP

#include "external-potential-impl.hpp"


namespace simploce {

    /**
     * Holds a collection of virtual planes. Particles inside a simulation box can
     * interact with these planes.
     * @see class VirtualPlane
     */
    class VirtualPlanes : public external_potential_impl {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @param spacing Spacing or distance between virtual planes. Default value is 2.0 nm.
         * @param eps_r Relative permittivity. Default is the value for water.
         */
        VirtualPlanes(box_ptr_t box,
                      bc_ptr_t bc,
                      dist_t spacing = 2.0,
                      real_t eps_r = 78.5);

        std::pair<energy_t, force_t>
        operator () (const p_ptr_t& particle) const override;

        /**
         * Sets the initial values of surface charge densities according to
         * the current positions of charged particles in the given particle system.
         */
        void
        initialize(const p_system_ptr_t& particleSystem) override;

        /**
         * Updates the values of surface charge densities according to
         * all current positions of charged particles in the given particle system.
         * Frozen particles are excluded.
         */
        void
        update(const p_system_ptr_t& particleSystem) override;

        /**
         * Updates values of surface charge densities according to
         * the current position of the gien particle.
         * @param particle Particle.
         */
        void
        update(const p_ptr_t& particle) override;

        /**
         * Restores the previous values of surface charge densities.
         */
        void
        fallback() override;

        /**
         * Returns virtual planes.
         * @return Virtual planes.
         */
        static std::vector<vplane_ptr_t>
        virtualPlanes() ;

        /**
         * Returns the total joint surface charge density, that is the sum of surface charge
         * densities over all virtual planes.
         * @return Total joint surface charge density.
         */
        static srf_charge_density_t
        surfaceChargeDensity() ;

        /**
         * Writes the location, averaged surface charge density, total accumulated charge
         * for each virtual plane to the output stream.
         * The name of the output file is 'vplanes.dat' in the 'current' directory.
         */
        void
        complete() const override;

    private:

        /**
         * Determines changes in state according to locations of particles in the given particle
         * system. Excludes frozen particles.
         */
        void
        determineStateChanges(const p_system_ptr_t& particleSystem);

        /**
         * Determines changes in state according to the given particle's location. Excludes frozen particles.
         */
        void
        determineStateChanges(const p_ptr_t& particle);

        /**
         * Determines changes in state according to the given particles' locations. Excludes frozen particles.
         */
        std::vector<charge_t>
        determineStateChanges(const std::vector<p_ptr_t>& particles);

        /**
         * Assigns surface charge densities to virtual planes.
         * @param updateAccumulated Update accumulated charge before resetting surface densities.
         */
        void resetSurfaceChargeDensities();

        box_ptr_t box_;
        bc_ptr_t bc_;
        dist_t spacing_;
        real_t eps_r_;
    };

}

#endif //SIMULATION_VIRTUAL_PLANES_HPP
