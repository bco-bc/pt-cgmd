/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/16/21.
 */

#include "simploce/potentials/lj.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/particle.hpp"

namespace simploce {

    LJ::LJ(ff_ptr_t forceField, bc_ptr_t bc) :
        forceField_(std::move(forceField)), bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    LJ::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        // Get C12 and C6 parameters.
        auto params = forceField_->lennardJones(p1->spec(), p2->spec());
        auto C12 = params.first;
        auto C6 = params.second;

        // Current positions.
        const auto &r1 = p1->position();
        const auto &r2 = p2->position();

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(r1, r2);
        auto Rij = norm<real_t>(rij);
        real_t Rij2 = Rij * Rij;

        // Forces and energy.
        return std::move(this->forceAndEnergy(rij, Rij, Rij2, C12, C6));
    }

    std::pair<energy_t, force_t>
    LJ::forceAndEnergy(const dist_vect_t& rij,
                       real_t Rij,
                       real_t Rij2,
                       real_t C12,
                       real_t C6) {
        // Potential/interaction energy, kJ/mol.
        real_t Rij6 = Rij2 * Rij2 * Rij2;
        real_t Rij12 = Rij6 * Rij6;
        real_t t1 = C12 / Rij12;
        real_t t2 = C6 / Rij6;
        energy_t energy{t1 - t2};

        // Force on particle #1, kJ/(mol nm).
        dist_vect_t unitVector = rij / Rij;
        real_t dLJdR = -6.0 * ( 2.0 * t1 - t2 ) / Rij;
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = -dLJdR * unitVector[k];
        }

        // Done.
        return std::move(std::make_pair(energy, f));

    }

}

