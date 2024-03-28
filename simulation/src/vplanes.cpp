/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 2 February 2024. Adapted from original version (1996).
 * @see https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(199612)17:16%3C1783::AID-JCC1%3E3.0.CO;2-J
 */

#include "simploce/potentials/vplanes.hpp"
#include "simploce/potentials/vplane.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/file.hpp"
#include <stdexcept>

namespace simploce {
    namespace virtual_planes {

        // All virtual planes.
        static std::vector<vplane_ptr_t> virtualPlanes_{};

        // Current charge values for each plane.
        static std::vector<charge_t> state_{};

        // Changes in state relative to the previous state.
        static std::vector<charge_t> difference_{};

        // Accumulated charge values for each plane.
        static std::vector<charge_t> accumulated_{};

        // Keeps track how often an update tool place.
        static int counter_{0};

        /**
         * Returns area of the xy-plane.
         * @param box Box
         * @return Area.
         */
        static area_t
        area(const box_ptr_t& box) {
            static area_t area = box->lengthX() * box->lengthY();
            return area;
        }

        /**
         * Increments counter by one.
         * @return Counter value after incrementCounter.
         */
        static int
        incrementCounter() {
            counter_ += 1;
            return counter_;
        }

        /**
         * Returns counter value.
         * @return Counter value.
         */
        static int
        counter() {
            return counter_;
        }

        /**
         * Updates state and accumulated.
         */
        static void
        updateStateAndAccumulated() {
            for (auto k = 0; k != accumulated_.size(); ++k) {
                state_[k] += difference_[k];
                accumulated_[k] += state_[k];
            }
        }

        /**
         * Removes current state from accumulated, restores previous state, and adds restored
         * state to accumulated.
         */
        static void
        revertStateAndAccumulated() {
            for (auto k = 0; k != state_.size(); ++k) {
                accumulated_[k] -= state_[k];  // Remove current state.
                state_[k] -= difference_[k];   // Restore previous state to become current state.
                accumulated_[k] += state_[k];  // Add previous state.
            }
        }

        /**
         * Sets up virtual planes.
         * @param box Box
         * @param bc Boundary conditions.
         * @param spacing Requested spacing between planes.
         * @param eps_r Relative permittivity.
         * @return Actual spacing between planes. May be slight different from requested spacing.
         */
        static dist_t
        setup(const box_ptr_t& box,
              const bc_ptr_t& bc,
              const dist_t& spacing,
              real_t eps_r) {
            auto numberOfPlanes = std::size_t(box->lengthZ() / spacing());
            auto distance = box->lengthZ() / real_t(numberOfPlanes);
            virtual_planes::virtualPlanes_ = std::vector<vplane_ptr_t>(numberOfPlanes);
            for (auto i = 0; i != numberOfPlanes; ++i) {
                auto location = i * distance;
                virtual_planes::virtualPlanes_[i] =
                        std::make_shared<VirtualPlane>(box, bc, location, eps_r);
            }
            virtual_planes::state_ = std::vector<charge_t>(numberOfPlanes, 0.0);
            virtual_planes::difference_ = std::vector<charge_t>(numberOfPlanes, 0.0);
            virtual_planes::accumulated_ = std::vector<charge_t>(numberOfPlanes, 0.0);
            return distance;
        }
    }

    VirtualPlanes::VirtualPlanes(box_ptr_t box,
                                 bc_ptr_t bc,
                                 dist_t spacing,
                                 real_t eps_r) :
    external_potential_impl{}, box_{std::move(box)}, bc_{std::move(bc)},
    spacing_{spacing}, eps_r_{eps_r} {
        util::Logger logger{"simploce::VirtualPlanes::VirtualPlanes()"};
        logger.trace("Entering.");

        // Validate.
        if (!box_) {
            throw std::domain_error("Box must be provided.");
        }
        if (!bc_) {
            throw std::domain_error("Boundary conditions must be provided.");
        }
        if (spacing_() == 0) {
            throw std::domain_error("There must be some space between virtual planes.");
        }
        if (eps_r_ <= 0) {
            throw std::domain_error("Relative permittivity must be a positive number.");
        }

        // Set up the virtual planes.
        spacing_ = virtual_planes::setup(box_, bc_, spacing_, eps_r_);

        // Log some information.
        auto numberOfPlanes = virtual_planes::virtualPlanes_.size();
        logger.info(util::to_string(spacing_) + ": Spacing between virtual planes.");
        logger.info(std::to_string(numberOfPlanes) + ": Number of virtual planes.");
        logger.info(std::to_string(eps_r_) + ": Relative permittivity.");

        logger.trace("Leaving.");
    }

    std::pair<energy_t, force_t>
    VirtualPlanes::operator()(const simploce::p_ptr_t &particle) const {
        util::Logger logger("simploce::VirtualPlanes::operator () ()");
        logger.trace("Entering.");

        logger.debug("This external potential is included.");
        energy_t energy{0.0};
        force_t force{0.0, 0.0, 0.0};
        auto& virtualPlanes = virtual_planes::virtualPlanes_;
        for  (const auto& plane: virtualPlanes) {
            auto result = plane->operator()(particle);
            energy += result.first;
            force += result.second;
        }
        logger.debug(util::to_string(energy) + ": Total interaction energy.");
        logger.debug(util::to_string(force) + ": Total force.");

        logger.trace("Leaving.");
        return std::move(std::make_pair(energy, force));
    }

    void
    VirtualPlanes::initialize(const p_system_ptr_t &particleSystem) {
        static util::Logger logger{"simploce::VirtualPlanes::initialize()"};
        logger.trace("Entering.");

        this->determineStateChanges(particleSystem);
        virtual_planes::updateStateAndAccumulated();
        virtual_planes::incrementCounter();
        this->resetSurfaceChargeDensities();
        logger.info(util::to_string(VirtualPlanes::surfaceChargeDensity()) +
        ": Joint total surface charge density of all virtual planes.");

        logger.trace("Leaving.");
    }

    void
    VirtualPlanes::update(const p_system_ptr_t &particleSystem) {
        static util::Logger logger("simploce::VirtualPlanes::update(const p_system_ptr_t &particleSystem)");
        logger.trace("Entering.");

        // Update.
        this->determineStateChanges(particleSystem);
        virtual_planes::updateStateAndAccumulated();
        virtual_planes::incrementCounter();
        this->resetSurfaceChargeDensities();

        logger.trace("Leaving");
    }

    void
    VirtualPlanes::update(const simploce::p_ptr_t &particle) {
        static util::Logger logger("simploce::VirtualPlanes::update(const p_ptr_t &particle)");
        logger.trace("Entering.");

        this->determineStateChanges(particle);
        virtual_planes::updateStateAndAccumulated();
        virtual_planes::incrementCounter();
        this->resetSurfaceChargeDensities();

        logger.trace("Leaving");
    }

    void
    VirtualPlanes::fallback() {
        // Restores previous state.
        util::Logger logger("simploce::VirtualPlanes::fallback()");
        logger.trace("Entering.");

        // No change of counter.
        virtual_planes::revertStateAndAccumulated();
        this->resetSurfaceChargeDensities();

        logger.trace("Leaving.");
    }

    std::vector<vplane_ptr_t>
    VirtualPlanes::virtualPlanes() {
        return virtual_planes::virtualPlanes_;
    }

    srf_charge_density_t
    VirtualPlanes::surfaceChargeDensity() {
        srf_charge_density_t total{0.0};
        for (const auto &p: virtual_planes::virtualPlanes_) {
            total += p->surfaceChargeDensity();
        }
        return total;
    }

    void
    VirtualPlanes::complete() const {
        util::Logger logger("simploce::VirtualPlanes::complete()");
        logger.trace("Entering.");

        std::ofstream stream;
        std::string fn = "vplanes.dat";
        util::open_output_file(stream, fn);
        for (auto k = 0; k != virtual_planes::virtualPlanes_.size(); ++k) {
            auto& p = virtual_planes::virtualPlanes_[k];
            stream << *p << virtual_planes::accumulated_[k] << std::endl;
        }
        stream.close();

        auto counter = virtual_planes::counter();
        logger.debug(std::to_string(counter) +
                    ": Number of times surface charge densities of virtual planes was reset.");
        logger.info(fn +
                    ": Virtual planes data (location, surface charge density) "
                    "was written to this output file.");
        logger.info(util::to_string(VirtualPlanes::surfaceChargeDensity()) +
                    ": Joint total surface charge density of all virtual planes.");

        logger.trace("Leaving.");
    }

    void
    VirtualPlanes::determineStateChanges(const p_system_ptr_t& particleSystem) {
        // Use all charge particles, but exclude frozen particles, to calculate the change
        // in the surface charge distribution over virtual planes.
        virtual_planes::difference_ =
                particleSystem->doWithAll<std::vector<charge_t>>([this] (
                        std::vector<p_ptr_t>& all) {
            auto difference = this->determineStateChanges(all);
            return std::move(difference);
        });
    }

    void
    VirtualPlanes::determineStateChanges(const simploce::p_ptr_t &particle) {
        std::vector<p_ptr_t> particles{particle};
        virtual_planes::difference_ = this->determineStateChanges(particles);
    }

    std::vector<charge_t>
    VirtualPlanes::determineStateChanges(const std::vector<p_ptr_t>& particles) {
        auto difference = std::vector<charge_t>(virtual_planes::virtualPlanes_.size(), 0.0);
        auto counter = virtual_planes::counter();
        for (const auto& p : particles) {
            auto Q = p->charge()();
            if (!p->frozen() && std::fabs(Q) > 0.0) {
                if (counter > 0) {
                    // Remove charge from virtual plane.
                    auto r_p = this->bc_->placeInside(p->previousPosition());
                    auto z_p = r_p[2];
                    auto index_p = std::size_t(z_p / this->spacing_());
                    difference[index_p] -= p->charge();
                }
                // Add charge to virtual plane.
                auto r = this->bc_->placeInside(p->position());
                auto z = r[2];
                auto index = std::size_t(z / this->spacing_());
                difference[index] += p->charge();
            }
        }
        return std::move(difference);
    }

    void VirtualPlanes::resetSurfaceChargeDensities() {
        static util::Logger logger("simploce::VirtualPlanes::resetSurfaceChargeDensities()");
        logger.trace("Entering.");

        auto& virtualPlanes = virtual_planes::virtualPlanes_;
        auto& accumulated = virtual_planes::accumulated_;
        auto area = virtual_planes::area(box_);
        auto counter = virtual_planes::counter() == 0 ? 1 : virtual_planes::counter();
        for (auto k = 0; k != virtualPlanes.size(); ++k) {
            auto average = accumulated[k]() / counter;      // Average charge.
            srf_charge_density_t sigma = average / area();
            auto& p = virtualPlanes[k];
            p->reset(sigma);
        }

        logger.trace("Leaving.");
    }

}
