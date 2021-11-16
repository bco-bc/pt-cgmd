/**
 * Author: André H. Juffer.
 * Created on 15/11/2021, 14:20
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_HALVE_ATTRACTIVE_QP_HPP
#define SIMULATION_HALVE_ATTRACTIVE_QP_HPP

#include "pair-potential.hpp"
#include "force-field.hpp"
#include "bc.hpp"
#include "simploce/particle/particle.hpp"
#include <utility>

namespace simploce {

    /**
     * Halve attractive quartic potential, U(r) = 0.5 * k * (r - r0)^4, where r is a distance,
     * k is a force constant (fc) and r0 is the equilibrium distance. Halve attractive implies
     * U(r) = 0 for r <= r0.
     * @see <a href="https://dx.doi.org/10.1063/1.3553378">
     * Riniker and van Gunsteren, J. Chem. Phys., 134, 084110.2011.
     * </a>
     * @tparam P Particle type.
     */
    template <typename P>
    class HalveAttractiveQP : public pair_potential<P> {
    public:

        using p_ptr_t = typename pair_potential<P>::p_ptr_t;

        HalveAttractiveQP(ff_ptr_t forceField, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
    };

    template <typename P>
    HalveAttractiveQP<P>::HalveAttractiveQP(ff_ptr_t forceField, bc_ptr_t bc) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)} {
    }

    template <typename P>
    std::pair<energy_t, force_t>
    HalveAttractiveQP<P>::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        // Get r0 and fc.
        auto params = forceField_->halveAttractiveQuarticParameters(p1->spec(),p2->spec());
        auto r0 = params.first;
        auto fc = params.second;

        // Current positions.
        const auto &r1 = p1->position();
        const auto &r2 = p2->position();

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(r1, r2);
        auto Rij = norm<real_t>(rij);
        real_t dR = Rij - r0;

        if ( dR > 0.0 ) {
            real_t dR3 = dR * dR * dR;
            real_t dR4 = dR * dR3;
            energy_t energy{0.5 * fc * dR4};
            real_t derQP = 2.0 * fc * dR3;
            dist_vect_t unitVector = rij / Rij;
            force_t f{};
            for (std::size_t k = 0; k != 3; ++k) {
                f[k] = -derQP * unitVector[k];
            }
            return  std::move(std::make_pair(energy, f));
        } else {
            return  std::move(std::make_pair(energy_t{0.0}, force_t{0.0, 0.0, 0.0}));
        }



    }
}

#endif //SIMULATION_HALVE_ATTRACTIVE_QP_HPP
