/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 15:14 PM
 */

#ifndef PARTICLES_BASE_CONTENT_HANDLER_HPP
#define PARTICLES_BASE_CONTENT_HANDLER_HPP

#include "content-handler.hpp"

namespace simploce {

    /**
     * Convenience implementation. Simply logs everything it receives.
     */
    class BaseContentHandler : public content_handler {
    public:

        BaseContentHandler();

        ~BaseContentHandler() override;

        void start(const std::string &title) override;

        void end() override;

        void startAtomGroup(const std::string &name) override;

        void endAtomGroup() override;

        void startMolecule(const std::string &name) override;

        void endMolecule() override;

        void startAtom(const std::string &name) override;

        void atomCoordinates(const position_t& r) override;

        void atomCharge(const charge_t& charge) override;

        void endAtom() override;

        void index(int index) override;

        /**
         * Returns number of atoms received.
         * @return Number of atoms.
         */
        std::size_t numberOfAtoms() const;

        /**
         * Returns number of atom groups received.
         * @return Number of atom groups.
         */
        std::size_t numberOfAtomGroups() const;

        /**
         * Returns number of molecules received.
         * @return Number of molecules.
         */
        std::size_t numberOfMolecules() const;

        /**
         * Returns total charge.
         * @return Total charge.
         */
        charge_t totalCharge() const;

    private:

        std::size_t numberOfAtoms_;
        std::size_t numberOfAtomGroups_;
        std::size_t numberOfMolecules_;
        charge_t totalCharge_;
    };

}

#endif //PARTICLES_BASE_CONTENT_HANDLER_HPP
