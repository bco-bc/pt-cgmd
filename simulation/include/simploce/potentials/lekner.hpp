/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 2 February 2024
 */
#ifndef SIMULATION_LEKNER_HPP
#define SIMULATION_LEKNER_HPP

#include "pair-potential.hpp"
#include "../conf/s-conf.hpp"

namespace simploce {

    /**
     * Lekner summation technique for calculating electrostatic interactions where there
     * is periodicity in the x- and y-direction, but -not- in the z-direction. Suitable
     * for simulations at surfaces. Currently only for Monte Carlo simulations as the
     * force calculation is not implemented.
     *
     * References:
     * 1. J. Lekner, "Summation of dipolar fields in simulated liquid-vapour interfaces.",
     *    Physica A, 157, 826-838, 1989.
     * 2. J. Lekner, "Summation of Coulomb fields in computer-simulated disordered systems",
     *    Physica A, 176, 485-498, 1991.
     * 3. Clark, A.T., Madden, T.J. and Warren, P.B., "Summation of electrostatic interactions
     *    in quasi-two-dimensional simulations.", Molecular Physics, 87, 1063-1069, 1996.
     * 4. N. Gronbech-Jensen, G. Hummer and Keith M. Beardmore, "Lekner summation of Coulomb
     *    interactions in partially periodic systems.", Molecular Physics., 92, 941-945, 1997.
     * 5. A.H. Juffer, C.M. Shepherd, H.J. Vogel, "Protein–membrane electrostatic interactions:
     *    Application of the Lekner summation technique.", J. Chem. Phys., 114, 1892 - 1905,
     *    2001. http://dx.doi.org/10.1063/1.1334901
     */
    class Lekner : public pair_potential {
    public:

        /**
         * Constructor
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @param eps Tolerance for convergence of sums in Eqs (7,8) of ref. 5.
         * @param n_max Maximum number of terms for the summation over n
         * (Eqs. (7,8) reference 5).
         * @param k_max Maximum number of terms for the summation over k
         * (Eqs. (7,8) reference 5).
         */
        explicit Lekner(box_ptr_t box,
                        bc_ptr_t bc,
                        real_t eps = conf::SMALL,
                        std::size_t n_max = 100,
                        std::size_t k_max = 100);

        /*
         * @return Energy (forces are always zero).
         */
        std::pair<energy_t, force_t> operator () (const p_ptr_t &pi,
                                                  const p_ptr_t &pj) override;

        /**
         * Interaction energy of particles.
         * @param Rij Distance vector.
         * @param qi Charge of particle #1.
         * @param qj Charge of particle #2.
         * @return Energy (forces are always zero).
         */
        std::pair<energy_t, force_t>
        forceAndEnergy(const dist_vect_t& Rij,
                       const charge_t& qi,
                       const charge_t& qj) const;

    private:

        box_ptr_t box_;
        bc_ptr_t bc_;
        real_t eps_;
        size_t n_max_;
        size_t k_max_;
    };

}

#endif //SIMULATION_LEKNER_HPP
