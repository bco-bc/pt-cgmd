/*
 * File:   protonation-site.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 4:43 PM
 */

#ifndef PROTONATION_SITE_HPP
#define PROTONATION_SITE_HPP

#include "particle-group.hpp"
#include "protonatable.hpp"
#include "particle.hpp"
#include "particle-spec.hpp"
#include "simploce/util/logger.hpp"
#include <stdexcept>
#include <utility>
#include <vector>
#include <memory>
#include <cassert>

namespace simploce {
    
    /**
     * A protonatable particle group. (De)protonation results in a redistribution of
     * charge and mass over the constituting particles. A maximum of one proton
     * can be bound or released.
     * @param T Particle type.
     */
    template <typename P>
    class ProtonationSite: public ParticleGroup {
    public:

        /**
         * Constructor. Creates a protonation site in a deprotonated state.
         * @param name Name, e.g. COOH.
         * @param particles Constituting particles. Must contain at least two particle.
         * @param bonds Holds identifiers of particles forming bonds. Particles
         * forming bonds must be constituting particles.
         * @param deprotonatedSpecs Specifications for the deprotonatedSpecs state. The
         * specification names must match those of the given constituting particles.
         * @param protonatedSpecs Specifications for the protonatedSpecs state. The
         * specification names must must match those of the given constituting particles.
         */
        ProtonationSite(std::string name,
                        const std::vector<p_ptr_t>& particles,
                        const std::vector<id_pair_t>& bonds,
                        const std::vector<spec_ptr_t>& deprotonatedSpecs,
                        const std::vector<spec_ptr_t>& protonatedSpecs);

        virtual ~ProtonationSite();
        
        /**
         * Returns site name.
         * @return Name.
         */
        std::string name() const;

        void protonate();

        void deprotonate();

        bool isProtonated() const;

        std::size_t protonationState() const;

    private:

        static util::Logger logger_;

        void matchSpecs_(const std::vector<spec_ptr_t>& deprotonatedSpecs,
                         const std::vector<spec_ptr_t>& protonatedSpecs);
        
        std::string name_;
        std::vector<spec_ptr_t> deprotonatedSpecs_;
        std::vector<spec_ptr_t> protonatedSpecs_;
        bool protonated_;
    };
    
    template <typename P>
    ProtonationSite<P>::ProtonationSite(std::string name,
                                        const std::vector<p_ptr_t>& particles,
                                        const std::vector<id_pair_t>& bonds,
                                        const std::vector<spec_ptr_t>& deprotonatedSpecs,
                                        const std::vector<spec_ptr_t>& protonatedSpecs) :
        ParticleGroup(particles, bonds), name_(std::move(name)),
        deprotonatedSpecs_{}, protonatedSpecs_{}, protonated_{true} {
        if ( name_.empty() ) {
            std::string message = "A name must be provided.";
            util::logAndThrow(ProtonationSite<P>::logger_, message);
        }
        std::size_t size = this->particles().size();
        if ( size <= 1 ) {
            std::string message = "Must consist of at least two particles.";
            util::logAndThrow(ProtonationSite<P>::logger_, message);
        }
        if (deprotonatedSpecs.size() != size || protonatedSpecs.size() != size ) {
            std::string message =
                    "Number of particles is not equal to number of particle specifications.";
            util::logAndThrow(ProtonationSite<P>::logger_, message);
        }
        this->matchSpecs_(deprotonatedSpecs, protonatedSpecs);
        this->deprotonate();
    }

    template <typename P>
    ProtonationSite<P>::~ProtonationSite<P>() = default;

    template <typename P>
    std::string
    ProtonationSite<P>::name() const {
        return name_;
    }
    
    template <typename P>
    void
    ProtonationSite<P>::protonate() {
        assert( !this->isProtonated() );
        auto particles = this->particles();
        for (auto particle : particles) {
            auto spec = particle->spec();
            bool found = false;
            for (auto protonatedSpec : protonatedSpecs_) {
                if (spec->name() == protonatedSpec->name() ) {
                    particle->resetSpec(protonatedSpec);
                }
            }
        }
        protonated_ = true;
    }
    
    template <typename P>
    void
    ProtonationSite<P>::deprotonate() {
        assert( this->isProtonated() );
        auto particles = ParticleGroup::particles();
        for (auto particle : particles) {
            auto spec = particle->spec();
            bool found = false;
            for (auto deprotonatedSpec : deprotonatedSpecs_) {
                if (spec->name() == deprotonatedSpec->name() ) {
                    particle->resetSpec(deprotonatedSpec);
                }
            }
        }
       protonated_ = false;
    }
    
    template <typename P>
    bool
    ProtonationSite<P>::isProtonated() const {
        return protonated_;
    }
    
    template <typename P>
    std::size_t
    ProtonationSite<P>::protonationState() const {
        return (this->isProtonated() ? 1 : 0);
    }

    template <typename P>
    util::Logger ProtonationSite<P>::logger_("ProtonationSite");
    
    template <typename P>
    void
    ProtonationSite<P>::matchSpecs_(const std::vector<spec_ptr_t>& deprotonatedSpecs,
                                    const std::vector<spec_ptr_t>& protonatedSpecs) {
        auto particles = this->particles();
        for (const auto& particle : particles) {
            auto name = particle->spec()->name();
            bool found = false;
            for (const auto& s: deprotonatedSpecs) {
                if (name == s->name()) {
                    this->deprotonatedSpecs_.push_back(s);
                    found = true;
                }
            }
            if ( !found ) {
                std::string message =
                        "Missing spec '" + name + "' for deprotonated state " +
                        "for site '" + name_ + "'.";
                util::logAndThrow(ProtonationSite<P>::logger_, message);
            }
            found = false;
            for (const auto& s: protonatedSpecs) {
                if (name == s->name()) {
                    this->protonatedSpecs_.push_back(s);
                    found = true;
                }
            }
            if ( !found ) {
                std::string message =
                        "Missing spec '" + name + "' for protonated state " +
                        "for site '" + name_ + "'.";
                util::logAndThrow(ProtonationSite<P>::logger_, message);
            }
        }
    }
     
}


#endif /* PROTONATION_SITE_HPP */

