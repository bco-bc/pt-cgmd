/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 29 February 2024.
 */

#ifndef SIMULATION_EXTERNAL_POTENTIAL_IMPL_HPP
#define SIMULATION_EXTERNAL_POTENTIAL_IMPL_HPP

#include "external-potential.hpp"

namespace simploce {

    /**
     * All external potential may be derived from this base class.
     */
    struct external_potential_impl : public external_potential {

        external_potential_impl() = default;

        ~external_potential_impl() override = default;

        /**
         * Does nothing.
         * @return Zero energy and force.
         */
        std::pair<energy_t, force_t>
        operator () (const p_ptr_t& particle) const override;

        /**
         * Does nothing.
         */
        void
        initialize(const p_system_ptr_t& particleSystem) override;

        /**
         * Does nothing.
         */
        void
        update(const p_system_ptr_t& particleSystem) override;

        /**
         * Does nothing.
         */
        void
        update(const p_ptr_t& particle) override;

        /**
         * Does nothing.
         */
        void
        fallback() override;

        /**
         * Does nothing.
         */
        void
        complete() const override;
    };
}

#endif //SIMULATION_EXTERNAL_POTENTIAL_IMPL_HPP
