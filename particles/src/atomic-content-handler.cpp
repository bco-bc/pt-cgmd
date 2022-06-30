/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on JUne 2, 2022.
 */

#include "simploce/particle/atomic-content-handler.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/util/logger.hpp"

#include <utility>

namespace simploce {

    struct AtomicContentHandlerHelper {
        util::Logger logger_{"simploce::AtomicContentHandler"};

        void start(const std::string &title);
        void startAtom(const std::string &name);
        void atomCoordinates(const position_t& r);
        void endAtom();
        void startAtomGroup(const std::string& name);
        void endAtomGroup();

        void particleSystem(const p_system_ptr_t& particleSystem);
        void catalog(const spec_catalog_ptr_t& catalog);
        void excludeWater(bool flagValue);

        p_system_ptr_t particleSystem_;
        spec_catalog_ptr_t catalog_;
        bool excludeWater_;

        bool isWater_{false};              // Current group is water.
        p_ptr_t atom_{};                   // Current atom.
        bool newAtomGroup_{false};         // Signals a atom group started.
        std::vector<p_ptr_t> atomGroup_{}; // Current atom group.
    };

    static AtomicContentHandlerHelper helper_{};

    AtomicContentHandler::AtomicContentHandler(p_system_ptr_t particleSystem,
                                               spec_catalog_ptr_t catalog,
                                               bool excludeWater) {
        helper_.particleSystem(particleSystem);
        helper_.catalog(catalog);
        helper_.excludeWater(excludeWater);
    }

    AtomicContentHandler::~AtomicContentHandler() = default;

    void
    AtomicContentHandler::start(const std::string &title) {
       helper_.start(title);
    }

    void
    AtomicContentHandler::end() {
    }

    void
    AtomicContentHandler::startAtomGroup(const std::string &name) {
        helper_.startAtomGroup(name);
    }

    void
    AtomicContentHandler::endAtomGroup() {
        helper_.endAtomGroup();
    }

    void
    AtomicContentHandler::startMolecule(const std::string &name) {
        this->startAtomGroup(name);
    }

    void
    AtomicContentHandler::endMolecule() {
        this->endAtomGroup();
    }

    void
    AtomicContentHandler::startAtom(const std::string &name) {
        helper_.startAtom(name);
    }

    void
    AtomicContentHandler::atomCoordinates(const position_t &r) {
        helper_.atomCoordinates(r);
    }

    void
    AtomicContentHandler::atomCharge(const charge_t &charge) {
    }

    void
    AtomicContentHandler::endAtom() {
        helper_.endAtom();
    }

    void
    AtomicContentHandler::index(int index) {
    }

    void
    AtomicContentHandlerHelper::start(const std::string &title) {
        logger_.info(title + ": New atomic particle system.");
    }

    void
    AtomicContentHandlerHelper::startAtomGroup(const std::string &name) {
        static std::vector<std::string> water = {"H2O", "HOH"};
        newAtomGroup_ = true;
        if (excludeWater_) {
            if (std::find(water.begin(), water.end(), name) == water.end()) {
                isWater_ = false;
            } else {
                isWater_ = true;
            }
        }
    }

    void
    AtomicContentHandlerHelper::endAtomGroup() {
        isWater_ = false;
        std::vector<id_pair_t> bonds{};
        if ( atomGroup_.size() > 1) {
            particleSystem_->addParticleGroup(atomGroup_, bonds);
        }
        newAtomGroup_ = false;
        atomGroup_.clear();
    }

    void
    AtomicContentHandlerHelper::startAtom(const std::string &name) {
        if (!isWater_ || !excludeWater_) {
            auto spec = catalog_->lookupByElementName(name);
            atom_ = particleSystem_->addParticle(name, spec);
            if (newAtomGroup_) {
                atomGroup_.emplace_back(atom_);
            }
        }
    }

    void AtomicContentHandlerHelper::atomCoordinates(const position_t &r) {
        if ( !isWater_ || !excludeWater_) {
            atom_->position(r);
        }
    }

    void
    AtomicContentHandlerHelper::endAtom() {
        atom_.reset();
    }

    void
    AtomicContentHandlerHelper::particleSystem(const p_system_ptr_t &particleSystem) {
        particleSystem_ = particleSystem;
    }

    void
    AtomicContentHandlerHelper::catalog(const spec_catalog_ptr_t &catalog) {
        catalog_ = catalog;
    }

    void AtomicContentHandlerHelper::excludeWater(bool flagValue) {
        excludeWater_ = flagValue;
    }

}