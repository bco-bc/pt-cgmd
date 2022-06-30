/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland
 *
 * Created on 10/22/21.
 */

#ifndef PARTICLES_PROTONATABLE_COARSE_GRAINED_HPP
#define PARTICLES_PROTONATABLE_COARSE_GRAINED_HPP

#include "coarse-grained.hpp"
#include "protonatable-bead.hpp"
#include "particle-spec.hpp"
#include "particle-spec-catalog.hpp"
#include "particle-group.hpp"
#include "p-factory.hpp"
#include "p-util.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/box.hpp"
#include "simploce/particle/p-properties.hpp"
#include <boost/algorithm/string.hpp>
#include <vector>
#include <memory>
#include <string>
#include <utility>

namespace simploce {

    /**
     * A particle system composed of both protonatable and regular (non-protonatable) beads.
     * @tparam S Protonation state type.
     */
    template <typename S>
    class ProtonatableCoarseGrained : public CoarseGrained {
    public:

        /**
         * Protonatable bead pointer type.
         */
        using prot_bead_ptr_t = std::shared_ptr<ProtonatableBead<S>>;

        /**
         * Protonatable coarse grained particle system pointer type.
         */
        using prot_cg_mod_ptr_t = std::shared_ptr<ProtonatableCoarseGrained<S>>;

        /**
         * Reads a protonatable coarse grained particle model from an input stream.
         * @param stream Input stream.
         * @param catalog Particle specifications catalog.
         * @return Protonatable coarse grained particle model.
         * @param stream
         * @param catalog
         * @return Protonatable coarse grained particle model.
         */
        static prot_cg_mod_ptr_t obtainFrom(std::istream& stream,
                                            const spec_catalog_ptr_t& catalog);

        ProtonatableCoarseGrained();

        // Noncopyable.
        ProtonatableCoarseGrained(const ProtonatableCoarseGrained&) = delete;
        ProtonatableCoarseGrained& operator = (const ProtonatableCoarseGrained&) = delete;

        // Movable.
        ProtonatableCoarseGrained(ProtonatableCoarseGrained&& pcg) noexcept;
        ProtonatableCoarseGrained& operator = (ProtonatableCoarseGrained&& pcg) noexcept;

        /**
         * Adds a protonatable bead to this particle system. All arguments are required.
         * @param S Protonation state type.
         * @param name Name (does not need to be unique).
         * @param spec Bead specification.
         * @return Protonatable bead added.
         */
        prot_bead_ptr_t addProtonatableBead(const std::string& name, const spec_ptr_t& spec);

        /**
         * Returns total number of bound protons.
         * @return Number, always >= 0.
         */
        std::size_t protonationState() const;

        /**
         * Perform a "task" with protonatable beads. The given task must expose the operator
         * <code>
         *  R operator () (std::vector<prot_bead_ptr_t>& protonatableBeads);
         * </code>
         * where the argument represents all protonatable beads of this particle system.
         * @tparam R Return type.
         * @tparam TASK TASK type.
         * @param task Task. This may be a lambda expression.
         * @return Result of type R. May be void.
         */
        template <typename R, typename TASK>
        R doWithProtonatableBeads(const TASK& task) { return task(protonatableBeads_); }

        /**
         * Returns number of protonatable beads.
         * @return Number.
         */
        int numberProtonatableBeads() const;

        bool isProtonatable() const override { return true; }

        /**
         * Replaces particle groups by a single particle of given specification. This method also assigns
         * new particle sequence indices to the particles.
         * @param spec Particle specification.
         * @param nReplace Number of group to be replaced.
         */
        void replaceGroupsByParticles(const spec_ptr_t& spec,
                                      int nReplace);

        /**
         * Returns protonatable bead with given identifier.
         * @param id Identifier.
         * @return Protonatable bead or nullptr if non-existent.
         */
        prot_bead_ptr_t findProtonatable(const id_t& id);

    private:

        p_ptr_t createParticle_(const id_t& id,
                                int index,
                                const std::string& name,
                                const spec_ptr_t& spec) override;

        std::vector<prot_bead_ptr_t> protonatableBeads_;
    };

    template <typename S>
    typename ProtonatableCoarseGrained<S>::prot_cg_mod_ptr_t
    ProtonatableCoarseGrained<S>::obtainFrom(std::istream& stream,
                                             const spec_catalog_ptr_t& catalog) {
        util::Logger logger{"simploce::ProtonatableCoarseGrained<S>::obtainFrom"};
        auto protonatableCoarseGrained = std::make_shared<ProtonatableCoarseGrained<S>>();
        protonatableCoarseGrained->parse(stream, catalog);
        logger.debug("Number of protonatable particles: " +
                      util::toString(protonatableCoarseGrained->numberProtonatableBeads()));
        return protonatableCoarseGrained;
    }

    template <typename S>
    ProtonatableCoarseGrained<S>::ProtonatableCoarseGrained() :
        CoarseGrained{}, protonatableBeads_{} {
    }

    template <typename S>
    ProtonatableCoarseGrained<S>::ProtonatableCoarseGrained(ProtonatableCoarseGrained&& pcg) noexcept {
        auto ptr = &pcg;
        auto& all = this->all();
        all.clear();
        auto p = ptr->all();
        all = std::move((&pcg)->all());
        auto& free = this->free();
        free = std::move((&pcg)->free());
        auto& groups = this->groups();
        groups = std::move((&pcg)->groups());
        box_ptr_t box = factory::box(pcg.box()->size());
        this->box(box);
        protonatableBeads_ = std::move(pcg.protonatableBeads_);
        (&pcg)->clear();
    }

    template <typename S>
    ProtonatableCoarseGrained<S>&
    ProtonatableCoarseGrained<S>::operator = (ProtonatableCoarseGrained&& pcg) noexcept {
        auto ptr = &pcg;
        auto& all = this->all();
        all.clear();
        auto p = ptr->all();
        all = std::move((&pcg)->all());
        auto& free = this->free();
        free = std::move((&pcg)->free());
        auto& groups = this->groups();
        groups = std::move((&pcg)->groups());
        box_ptr_t box = factory::box(pcg.box()->size());
        this->box(box);
        protonatableBeads_ = std::move(pcg.protonatableBeads_);
        (&pcg)->clear();
        return *this;
    }

    template <typename S>
    typename ProtonatableCoarseGrained<S>::prot_bead_ptr_t
    ProtonatableCoarseGrained<S>::addProtonatableBead(const std::string& name,
                                                      const spec_ptr_t &spec) {
        util::Logger logger{"simploce::ProtonatableCoarseGrained<S>::addProtonatableBead"};
        if ( !spec->isProtonatable() ) {
            util::logAndThrow(logger,
                              spec->name() +
                              ": Not a protonatable particle specification for particle '"  +
                              name + "'. Use method addBead(...) instead.");
        }
        auto bead = this->addParticle(name, spec);
        return this->findProtonatable(bead->id());
    }

    template <typename S>
    std::size_t
    ProtonatableCoarseGrained<S>::protonationState() const {
        std::size_t n = 0;
        for (const auto& p : protonatableBeads_) {
            if ( p->isProtonated() ) {
                n += 1;
            }
        }
        return n;
    }

    template <typename S>
    int
    ProtonatableCoarseGrained<S>::numberProtonatableBeads() const {
        return protonatableBeads_.size();
    }

    template <typename S>
    void
    ProtonatableCoarseGrained<S>::replaceGroupsByParticles(const spec_ptr_t& spec,
                                                           int nReplace) {
        for (int i = 0; i != nReplace; ++i) {
            auto v = util::random();
            int index = int(v * this->groups().size());
            auto group = this->groups()[index];
            position_t r = group->position();
            this->removeGroup(group);
            for ( auto& p: group->particles() ) {
                this->remove(p);
            }
            std::string name = spec->name() + util::toString(this->numberOfParticles());
            if (spec->isProtonatable() ) {
                auto p = this->addProtonatableBead(name, spec);
                p->position(r);
            } else {
                auto p = this->addBead(name, spec);
                p->position(r);
            }
        }
        this->resetIndex();
    }

    template<typename S>
    typename ProtonatableCoarseGrained<S>::prot_bead_ptr_t
    ProtonatableCoarseGrained<S>::findProtonatable(const id_t& id) {
        return util::find(id, protonatableBeads_);
    }

    template<typename S>
    p_ptr_t
    ProtonatableCoarseGrained<S>::createParticle_(const id_t& id,
                                                  int index,
                                                  const std::string& name,
                                                  const spec_ptr_t& spec) {
        if ( spec->isProtonatable() ) {
            prot_bead_ptr_t bead = ProtonatableBead<S>::create(id, index, name, spec);
            protonatableBeads_.emplace_back(bead);
            return bead;
        } else {
            return Bead::create(id, index, name, spec);
        }
    }

    /**
     * Writes protonatable coarse grained particle model to output stream.
     * @param stream Output stream.
     * @param cg Coarse grained particle model.
     * @return Output stream.
     */
    template <typename S>
    std::ostream&
    operator << (std::ostream& stream, const ProtonatableCoarseGrained<S>& cg) {
        cg.write(stream);
        return stream;
    }

}

#endif //PARTICLES_PROTONATABLE_COARSE_GRAINED_HPP
