/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 13, 2019, 3:22 PM
 */

#ifndef FORCE_FIELD_HPP
#define FORCE_FIELD_HPP

#include "pair-lists.hpp"
#include "s-types.hpp"
#include "simploce/util/map2.hpp"
#include <string>
#include <utility>
#include <iostream>

namespace simploce {

    /**
     * The set of interaction parameters for computing forces on particles.
     */
    class ForceField {
    public:

        /**
         * Holds parameters for an interaction.
         */
        using int_spec_t = struct {
            std::string type;
            spec_ptr_t spec1;
            spec_ptr_t spec2;
            real_t C12;
            real_t C6;
            real_t fc;     // Force constant.
            real_t r0;     // Equilibrium distance.
        };

        /**
         * Reads an force field from an input stream.
         * @param stream Input stream.
         * @param catalog Particle specification catalog.
         * @return Force field.
         */
        static ff_ptr_t obtainFrom(std::istream& stream, const spec_catalog_ptr_t& catalog);

        ForceField();
        
        virtual ~ForceField();

        // Noncopyable
        ForceField(const ForceField&) = delete;
        ForceField& operator = (const ForceField&) = delete;

        // Movable.
        ForceField(ForceField&& forceField) noexcept ;
        ForceField& operator = (ForceField&& forceField) noexcept ;

        /**
         * Adds interaction specification to force field.
         * @param spec Interaction specification.
         */
        void addInteractionSpecification(const int_spec_t &spec);

        /**
         * Returns interaction parameter set.
         * @return Interaction parameters.
         */
        const std::vector<int_spec_t> &interactionSpecifications() const;

        /**
         * Returns Lennard-Jones interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return C12, C6 respectively.
         */
        std::pair<real_t, real_t> lennardJonesParameters(const spec_ptr_t &spec1,
                                                         const spec_ptr_t &spec2) const;

        /**
         * Returns Lennard-Jones plus Reaction Field electrostatic interaction parameters.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return C12, C6, eps_inside_rc, eps_outside_rc, respectively. (eps = relative permittivity,
         * rc = cutoff distance.)
         */
        std::tuple<real_t, real_t, real_t, real_t>
        lennardJonesReactionFieldParameters(const spec_ptr_t &spec1,
                                           const spec_ptr_t &spec2) const;

        /**
         * Returns harmonic bond interaction parameters for the given pair of
         * particle specifications
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return Equilibrium distance (r0) and force constant (fc).
         */
        std::pair<real_t, real_t> harmonicParameters(const spec_ptr_t &spec1,
                                                     const spec_ptr_t &spec2) const;

        /**
         * Returns halve attractive quartic bond interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return Equilibrium distance (r0) and force constant (fc).
         */
        std::pair<real_t, real_t> halveAttractiveQuarticParameters(const spec_ptr_t &spec1,
                                                                   const spec_ptr_t &spec2) const;

        /**
         * Sets value for the relative permittivity inside cutoff distance.
         * @param eps Relative permittivity.
         */
        void relativePermittivityInsideCutoff(real_t eps);

        /**
         * Returns value of inside relative permittivity. Its default value is 2.5.
         * @return Relative permittivity.
         */
        real_t relativePermittivityInsideCutoff() const;

        /**
         * Relative permittivity outside (beyond) cutoff distance. Its default value
         * is 78.5.
         * @return Relative permittivity.
         */
        real_t relativePermittivityBeyondCutoff() const;

        /**
         * Sets relative permittivity outside (beyond) cutoff distance.
         * @param eps Relative permittivity.
         */
        void relativePermittivityBeyondCutoff(real_t eps);

    private:

        real_t eps_inside_rc_;
        real_t eps_beyond_rc_;
        std::vector<int_spec_t> interactionsSpecs_;

    };

    /**
     * Write force field to an output stream.
     * @param stream Output stream.
     * @param forceField Force field.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const ForceField& forceField);
}

#endif /* FORCE_FIELD_HPP */

