/*
 * File:   properties.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 3:42 PM
 */

#ifndef P_PROPERTIES_HPP
#define P_PROPERTIES_HPP

#include "p-types.hpp"
#include "atom.hpp"
#include "bead.hpp"
#include <vector>

namespace simploce {
    /**
     * Properties of a collection of particles.
     */
    namespace properties {
        
        /**
         * Finds particle in a given particle collection.
         * @param id Particle identifier.
         * @param particles Particles.
         * @return Particle or nullptr if the particle cannot be identified.
         */
        template <typename P, template<typename, typename ...> class CONT = std::vector>
        std::shared_ptr<P> find(simploce::id_t id, const CONT<std::shared_ptr<P>>& particles)
        {
            for (auto p : particles) {
                if ( p->id() == id ) {
                    return p;
                }
            }
            return nullptr;
        }                    
        
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
    }
}

#endif /* P_PROPERTIES_HPP */

