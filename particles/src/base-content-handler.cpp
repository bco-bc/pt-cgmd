/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 15:14 PM
 */

#include "simploce/chem/base-content-handler.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"

namespace simploce {

    static util::Logger logger_{"simploce::BaseContentHandler"};

    BaseContentHandler::BaseContentHandler() :
            numberOfAtoms_{0}, numberOfAtomGroups_{0}, numberOfMolecules_{0}, totalCharge_{0.0} {
    }

    BaseContentHandler::~BaseContentHandler() = default;

    void
    BaseContentHandler::start(const std::string &title) {
        logger_.info("Start of chemical content: " + title);
    }

    void
    BaseContentHandler::end() {
        logger_.info("End of chemical content.");
        logger_.info("Number of atoms: " + util::to_string(this->numberOfAtoms()));
        logger_.info("Number of atom groups: " + util::to_string(this->numberOfAtomGroups()));
        logger_.info("Number of molecules: " + util::to_string(this->numberOfMolecules()));
        logger_.info("Total charge: " + util::to_string(this->totalCharge()));
    }

    void
    BaseContentHandler::startAtomGroup(const std::string &name) {
        logger_.info("Start of atom group: " + name);
        numberOfAtomGroups_ += 1;
    }

    void
    BaseContentHandler::endAtomGroup() {
        logger_.info("End of atom group.");
    }

    void
    BaseContentHandler::startMolecule(const std::string &name) {
        logger_.info("Start of molecule: " + name);
        numberOfMolecules_ += 1;
    }

    void
    BaseContentHandler::endMolecule() {
        logger_.info("End of molecule.");
    }

    void
    BaseContentHandler::startAtom(const std::string &name) {
        logger_.info("Start of atom: " + name);
        numberOfAtoms_ += 1;
    }

    void
    BaseContentHandler::atomCoordinates(const position_t& r) {
        logger_.info("Atom coordinates: " + util::to_string(r));
    }

    void
    BaseContentHandler::atomCharge(const charge_t& charge) {
        logger_.info("Atom charge: " + util::to_string(charge));
        totalCharge_ += charge;
    }

    void
    BaseContentHandler::endAtom() {
        logger_.info("End of atom.");
    }

    std::size_t
    BaseContentHandler::numberOfAtoms() const {
        return numberOfAtoms_;
    }
    
    void
    BaseContentHandler::index(int index) {
        logger_.info("Index received: " + util::to_string(index));
    }

    std::size_t
    BaseContentHandler::numberOfAtomGroups() const {
        return numberOfAtomGroups_;
    }

    std::size_t
    BaseContentHandler::numberOfMolecules() const {
        return numberOfMolecules_;
    }

    charge_t
    BaseContentHandler::totalCharge() const {
        return totalCharge_;
    }

}