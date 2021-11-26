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
            std::string typeName;   // Type name.
            spec_ptr_t spec1;       // First particle specification.
            spec_ptr_t spec2;       // Second particle specification.
            real_t C12;             // C12 Lennard-Jones interaction.
            real_t C6;              // C6 Lennard-Jones interaction.
            real_t fc;              // Force constant.
            real_t r0;              // Equilibrium distance.
            real_t eps_inside_rc;   // Relative permittivity inside cutoff distance.
            real_t eps_outside_rc;  // Relative permittivity outside cutoff distance.
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
         * Adds a non-bonded interaction specification to force field.
         * @param spec Interaction specification.
         */
        void addNonBondedSpecification(const int_spec_t& spec);

        /**
         * Adds a bonded interaction specification to force field.
         * @param spec Interaction specification.
         */
        void addBondedSpecification(const int_spec_t& spec);

        /**
         * Returns non-bonded interaction specifications.
         * @return Specifications.
         */
        const std::vector<int_spec_t>& nonBondedSpecifications() const;

        /**
         * Returns bonded interaction specifications.
         * @return Specifications.
         */
        const std::vector<int_spec_t>& bondedSpecifications() const;

        /**
         * Returns non-bonded Lennard-Jones interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return C12, C6 respectively.
         */
        std::pair<real_t, real_t> lennardJones(const spec_ptr_t &spec1,
                                               const spec_ptr_t &spec2) const;

        /**
         * Returns non-bonded reaction field electrostatic interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return eps_inside_rc, eps_outside_rc, respectively.
         * (eps = relative permittivity, rc = cutoff distance.)
         */
        std::pair<real_t, real_t> reactionField(const spec_ptr_t &spec1,
                                                const spec_ptr_t &spec2) const;


        /**
         * Returns non-bonded hard sphere reaction field electrostatic interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return eps_inside_rc, eps_outside_rc, respectively.
         * (eps = relative permittivity, rc = cutoff distance.)
         */
        std::pair<real_t, real_t> hardSphereReactionField(const spec_ptr_t &spec1,
                                                          const spec_ptr_t &spec2) const;
        /**
         * Returns non-bonded screened Coulomb electrostatic interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return eps_inside_rc (eps = relative permittivity, rc = cutoff distance.)
         */
        real_t screenedCoulomb(const spec_ptr_t &spec1,
                               const spec_ptr_t &spec2) const;

        /**
         * Returns hard-sphere screened Coulomb electrostatic interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return eps_inside_rc (eps = relative permittivity, rc = cutoff distance.)
         */
        real_t hardSphereScreenedCoulomb(const spec_ptr_t &spec1,
                                         const spec_ptr_t &spec2) const;

        /**
         * Returns shifted force electrostatic interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return eps_inside_rc (eps = relative permittivity, rc = cutoff distance.)
         */
        real_t shiftedForce(const spec_ptr_t &spec1,
                            const spec_ptr_t &spec2) const;

        /**
         * Returns hard-sphere shifted force electrostatic interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return eps_inside_rc (eps = relative permittivity, rc = cutoff distance.)
         */
        real_t hardSphereShiftedForce(const spec_ptr_t &spec1,
                                      const spec_ptr_t &spec2) const;


        /**
         * Returns non-bonded Lennard-Jones plus Reaction Field electrostatic interaction parameters.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return C12, C6, eps_inside_rc, eps_outside_rc, respectively.
         * (eps = relative permittivity, rc = cutoff distance.)
         */
        std::tuple<real_t, real_t, real_t, real_t>
        lennardJonesReactionField(const spec_ptr_t &spec1,
                                  const spec_ptr_t &spec2) const;

        /**
         * Returns bonded harmonic interaction parameters for the given pair of
         * particle specifications
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return Equilibrium distance (r0) and force constant (fc).
         */
        std::pair<real_t, real_t> harmonic(const spec_ptr_t &spec1,
                                           const spec_ptr_t &spec2) const;

        /**
         * Returns bonded halve attractive quartic interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return Equilibrium distance (r0) and force constant (fc).
         */
        std::pair<real_t, real_t> halveAttractiveQuartic(const spec_ptr_t &spec1,
                                                         const spec_ptr_t &spec2) const;

    private:

        std::vector<int_spec_t> nonBondedSpecs_;
        std::vector<int_spec_t> bondedSpecs_;

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

