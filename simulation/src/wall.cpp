/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 12/20/21.
 */

#include "simploce/potentials/wall.hpp"
#include "simploce/potentials/uniform-surface-charge-density.hpp"

#include <utility>
#include "simploce/potentials/lj.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/util/logger.hpp"

namespace simploce {

    Wall::Wall(real_t C12, real_t C6, bc_ptr_t bc, FlatSurface flatSurface, srf_charge_density_t sigma) :
            C12_{C12}, C6_{C6}, bc_{std::move(bc)}, flatSurface_{std::move(flatSurface)}, sigma_{sigma} {
    }

    std::pair<energy_t, force_t>
    Wall::operator () (const p_ptr_t& particle) const {
        static util::Logger logger{"simploce::Wall::operator () "};
        logger.debug("This external potential is included.");

        static UniformSurfaceChargeDensity uscg(sigma_, flatSurface_, 1.0, bc_, 0.0, false);

        // Interaction with charge surface.
        auto result1 = uscg.operator()(particle);

        // LJ interaction.
        const position_t r_out = particle->position();
        auto r = bc_->placeInside(r_out);
        auto pair = flatSurface_.distanceTo(r);
        auto Rij = pair.first;
        real_t Rij2 = Rij * Rij;
        auto rij = pair.second;
        auto result2 = LJ::forceAndEnergy(rij, Rij(), Rij2, this->C12_, this->C6_);

        return std::move(std::make_pair<energy_t, force_t>(result1.first + result2.first, result1.second + result2.second));
    }

}