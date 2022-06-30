/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on JUne 2, 2022.
 */

#ifndef PARTICLES_ATOMIC_CONTENT_HANDLER_HPP
#define PARTICLES_ATOMIC_CONTENT_HANDLER_HPP

#include "p-types.hpp"
#include "../chem/content-handler.hpp"

namespace simploce {

    /**
     * Creates an atomic particle system from chemical content. Water from the original
     * input source is excluded by default.
     */
    class AtomicContentHandler : public content_handler {
    public:

        /**
         * Constructor.
         * @param particleSystem Empty particle system.
         * @param catalog Particle specifications catalog.
         * @param excludeWater Exclude water? True by default.
         */
        AtomicContentHandler(p_system_ptr_t particleSystem,
                             spec_catalog_ptr_t catalog,
                             bool excludeWater = true);

        ~AtomicContentHandler() override;

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

    };
}

#endif //PARTICLES_ATOMIC_CONTENT_HANDLER_HPP
