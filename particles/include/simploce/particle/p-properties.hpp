/*
 * File:   properties.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 3:42 PM
 */

#ifndef P_PROPERTIES_HPP
#define P_PROPERTIES_HPP

#include "p-types.hpp"
#include "particle.hpp"
#include "atom.hpp"
#include "bead.hpp"
#include <vector>

namespace simploce {
    /**
     * Properties of a collection of particles.
     */
    namespace properties {
        
        /**
         * Returns total mass of a collection of particles.
         * @param P Particle type.
         * @param CONT Container type.
         * @param particles Particles.
         * @return Total mass.
         */
        template <typename P, template<typename, typename ...> class CONT = std::vector>
        mass_t mass(const CONT<std::shared_ptr<P>>& particles)
        {
            mass_t total{0.0};
            for (auto p : particles) {
                total += p->mass();
            }
            return total;
        }
           
        /**
         * Returns total charge of a collection of particles.
         * @param P Particle type.
         * @param CONT Container type.
         * @param particles Particles.
         * @return Total charge.
         */
        template <typename P, template<typename, typename ...> class CONT = std::vector>
        charge_t charge(const CONT<std::shared_ptr<P>>& particles)
        {
            charge_t total{0.0};
            for (auto p : particles) {
                total += p->charge();
            }
            return total;
        }
        
        /**
         * Returns the center of mass of a collection of particles.
         * @param P Particle type.
         * @param CONT Container type.
         * @param particles Particles.
         * @return Center of mass.
         */
        template <typename P, template <typename, typename...> class CONT = std::vector> 
        position_t centerOfMass(const CONT<std::shared_ptr<P>>& particles)
        {
            mass_t total{0.0};
            position_t cm{};
            for (auto p : particles) {
                total += p->mass();
                cm += p->mass()() * p->position();
            }
            return position_t{cm / total};
        }

        /**
         * Returns total linear momentum of a collection of particles.
         * @tparam P Particle type.
         * @tparam CONT Container type.
         * @param particles Particles.
         * @return Total linear momentum.
         */
        template <typename P, template <typename, typename...> class CONT = std::vector>
        momentum_t linearMomentum(const CONT<std::shared_ptr<P>>& particles) {
            momentum_t p{0.0, 0.0, 0.0};
            for (auto& particle : particles) {
                auto mass = particle->mass();
                auto v = particle->velocity();
                for (int k = 0; k != 3; ++k) {
                    p[k] += mass() * v[k];
                }
            }
            return p;
        }
    }
}

#endif /* P_PROPERTIES_HPP */

