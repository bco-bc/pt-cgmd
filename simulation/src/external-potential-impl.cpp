/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 29 February 2024.
 */

#include "simploce/potentials/external-potential-impl.hpp"

namespace simploce {

    std::pair<energy_t, force_t>
    external_potential_impl::operator()(const simploce::p_ptr_t &particle) const {
        return std::move(std::make_pair(energy_t{}, force_t{}));
    }

    void
    external_potential_impl::initialize(const p_system_ptr_t &particleSystem) {
    }

    void
    external_potential_impl::update(const p_system_ptr_t &particleSystem) {
    }

    void
    external_potential_impl::update(const simploce::p_ptr_t &particle) {
    }

    void
    external_potential_impl::fallback() {
    }

    void
    external_potential_impl::complete() const {
    }
}
