/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 13, 2019, 3:22 PM
 */

#ifndef FORCE_FIELD_HPP
#define FORCE_FIELD_HPP

#include "simploce/simulation/pair-list.hpp"
#include "simploce/types/s-types.hpp"
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
         * Holds parameters for a single interaction.
         *
         */
        using int_spec_t = struct {
            std::string typeName;   // Type name.
            spec_ptr_t spec1;       // First particle specification.
            spec_ptr_t spec2;       // Second particle specification.
            real_t C12;             // C12 Lennard-Jones interaction.
            real_t C6;              // C6 Lennard-Jones interaction.
            real_t fc;              // Force constant.
            real_t r0;              // Equilibrium distance.
            real_t eps_inside_rc;   // Relative permittivity inside cutoff distance for electrostatic interactions.
            real_t eps_outside_rc;  // Relative permittivity outside cutoff distance for electrostatic interactions.
            real_t eps_r;           // Relative permittivity.
            real_t deltaV;          // Electric potential difference.
            real_t distance;        // Distance between two points.
            el_field_t e0;          // External electric field applied in a direction.
            char direction;         // Direction (or axis) along which something is applied.
            std::string plane;      // Plane specification.
            real_t sigma;           // Uniform surface charge density.
            real_t delta;           // Stern layer for uniform surface density.
            real_t sigma1;          // Width Gaussian charge density.
            real_t sigma2;          // Width Gaussian charge density.
            real_t max_a;           // Maximum repulsion (MVV_DPD).
            force_t fe;             // External force applied in a direction.
            real_t spacing;         // Spacing between things.
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
         * Adds an external potential specification.
         * @param spec External potential specification.
         */
        void addExternalSpecification(const int_spec_t& spec);

        /**
         * Returns non-bonded interaction potential specifications.
         * @return Specifications.
         */
        const std::vector<int_spec_t>& nonBondedSpecifications() const;

        /**
         * Returns bonded interaction potential specifications.
         * @return Specifications.
         */
        const std::vector<int_spec_t>& bondedSpecifications() const;

        /**
         * Returns external potential specifications.
         * @return External potential specifications.
         */
        const std::vector<int_spec_t>& externalSpecifications() const;

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
         * @return eps_inside_rc, eps_outside_rc, respectively. (eps = relative permittivity, rc = cutoff distance.)
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
         * @see #reactionField(const spec_ptr_t &spec1, const spec_ptr_t &spec2)
         */
        std::pair<real_t, real_t> hardSphereReactionField(const spec_ptr_t &spec1,
                                                          const spec_ptr_t &spec2) const;
        /**
         * Returns non-bonded screened Coulomb (1/(eps * r)) electrostatic interaction parameters for the given pair of
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
         * @see #screenedCoulomb(const spec_ptr_t &spec1, const spec_ptr_t &spec2) const;
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
         * @return eps_inside_rc (eps = relative permittivity)
         * @see #shiftedForce(const spec_ptr_t &spec1, const spec_ptr_t &spec2) const
         */
        real_t hardSphereShiftedForce(const spec_ptr_t &spec1,
                                      const spec_ptr_t &spec2) const;

        /**
         * Returns interaction parameters for the electrostatic interaction between two overlapping reset densities.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return sigma1, sigma2, respectively. Sigma represents the "width" of each density.
         */
        std::pair<real_t, real_t>
        gaussianChargeDensity(const spec_ptr_t &spec1,
                              const spec_ptr_t &spec2) const;

        /**
         * Returns relative permittivity or dielectric constant for the hard-sphere 2D Lekner summation.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         */
        real_t
        hardSphereLekner(const spec_ptr_t &spec1,
                         const spec_ptr_t &spec2) const;

        /**
         * Returns interaction parameters for the combined interaction for
         * overlapping Gaussian reset densities and the soft repulsion. This is for
         * mesoscopic simulations.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return (sigma1, sigma2, max_a, respectively. Sigma represents the "width"
         * of each density, and max_a is the maximum strength of the soft repulsion
         * potential.
         */
        std::tuple<real_t, real_t, real_t>
        gaussianChargeDensitySoftRepulsion(const spec_ptr_t &spec1,
                                           const spec_ptr_t &spec2) const;


        /**
         * Returns non-bonded Lennard-Jones plus reaction field electrostatic interaction parameters.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return C12, C6, eps_inside_rc, eps_outside_rc, respectively.
         * (eps = relative permittivity, rc = cutoff distance.)
         */
        std::tuple<real_t, real_t, real_t, real_t>
        lennardJonesReactionField(const spec_ptr_t &spec1,
                                  const spec_ptr_t &spec2) const;

        /**
         * Returns non-bonded Lennard-Jones plus Shifted Force electrostatic interaction
         * parameters.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return C12, C6, eps_inside_rc, respectively.
         * (eps = relative permittivity, rc = cutoff distance.)
         */
        std::tuple<real_t, real_t, real_t>
        lennardJonesShiftedForce(const spec_ptr_t &spec1,
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
         * Returns  bonded halve attractive harmonic potential interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return Threshold distance (r0) and force constant (fc), respectively.
         */
        std::pair<real_t, real_t> halveAttractiveHarmonic(const spec_ptr_t &spec1,
                                                          const spec_ptr_t &spec2) const;

        /**
         * Returns bonded halve attractive quartic interaction parameters for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return Threshold distance (r0) and force constant (fc), respectively.
         */
        std::pair<real_t, real_t> halveAttractiveQuartic(const spec_ptr_t &spec1,
                                                         const spec_ptr_t &spec2) const;

        /**
         * Returns non-bonded interaction parameters for the conservative soft-repulsive potential for the given pair of
         * particle specifications.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return Maximum repulsion (a).
         */
        real_t softRepulsion(const spec_ptr_t &spec1,
                             const spec_ptr_t &spec2) const;

        /**
         * Returns bonded interaction potential parameters for a harmonic potential plus soft repulsive potential.
         * @param spec1 Particle specification #1.
         * @param spec2 Particle specification #2.
         * @return Equilibrium distance (r0), force constant (fc), Maximum repulsion (a).
         */
        std::tuple<real_t, real_t, real_t>
        harmonicSoftRepulsion(const spec_ptr_t &spec1,
                              const spec_ptr_t &spec2) const;

        /**
         * Returns external potential parameters for a uniform surface reset density.
         * @return Surface reset density, delta, eps_inside_rc (delta = Width Stern layer,
         * eps = relative permittivity, rc = cutoff distance), plane.
         */
        std::tuple<srf_charge_density_t, dist_t, real_t, std::string>
        uniformSurfaceChargeDensity() const;

    private:

        std::vector<int_spec_t> nonBondedSpecs_;
        std::vector<int_spec_t> bondedSpecs_;
        std::vector<int_spec_t> external_;

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

