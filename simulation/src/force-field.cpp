/*
 * Created by André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#include "simploce/simulation/force-field.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/s-types.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <boost/algorithm/string.hpp>
#include <utility>
#include <vector>
#include <memory>

namespace simploce {

    ff_ptr_t
    ForceField::obtainFrom(std::istream &stream,
                           const spec_catalog_ptr_t &catalog) {
        util::Logger logger("simploce::ForceField::obtainFrom");

        ff_ptr_t forceField = std::make_shared<ForceField>();
        char charBuffer[100];
        std::string stringBuffer;
        std::size_t numberOfInteractionTypes{0};

        // Skip two headers.
        std::getline(stream, stringBuffer);
        std::getline(stream, stringBuffer);

        // Read non-bonded interaction parameters.
        int nNonBonded;
        stream >> nNonBonded;
        std::getline(stream, stringBuffer);
        for (int k = 0; k != nNonBonded; ++k) {
            // Read interaction type.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string interactionType = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(interactionType);

            // Read particle specifications.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName1 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName1);
            auto spec1 = catalog->lookup(specName1);
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName2 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName2);
            auto spec2 = catalog->lookup(specName2);

            // Read interaction parameters.
            if ( interactionType == conf::LJ || interactionType == conf::LJ_RF ) {
                real_t C12, C6;
                stream >> C12 >> C6;
                int_spec_t spec;
                spec.type = interactionType;
                spec.spec1 = spec1;
                spec.spec2 = spec2;
                spec.C12 = C12;
                spec.C6 = C6;
                forceField->addInteractionSpecification(spec);
                numberOfInteractionTypes += 1;
            } else {
                logger.warn("Skipping interaction type '" + interactionType + "'.");
            }
            std::getline(stream, stringBuffer);
        }

        // Skip two headers.
        std::getline(stream, stringBuffer);
        std::getline(stream, stringBuffer);

        int nBonded;
        stream >> nBonded;
        std::getline(stream, stringBuffer);
        for (int k = 0; k != nBonded; ++k) {
            // Read interaction type.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string interactionType = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(interactionType);

            // Read particle specifications.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName1 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName1);
            auto spec1 = catalog->lookup(specName1);
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName2 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName2);
            auto spec2 = catalog->lookup(specName2);

            if ( interactionType == conf::HP || interactionType == conf::HA_QP ) {
                real_t r0, fc;
                stream >> r0 >> fc;
                int_spec_t spec;
                spec.type = interactionType;
                spec.spec1 = spec1;
                spec.spec2 = spec2;
                spec.r0 = r0;
                spec.fc = fc;
                forceField->addInteractionSpecification(spec);
                numberOfInteractionTypes += 1;
            } else {
                    logger.warn("Skipping interaction type '" + interactionType + "'.");
            }
            std::getline(stream, stringBuffer);
        }

        // Log some information.
        logger.debug("Number of interaction types: " + util::toString(numberOfInteractionTypes));

        return std::move(forceField);
    }

    ForceField::ForceField() : eps_inside_rc_{2.5}, eps_beyond_rc_{78.5}, interactionsSpecs_{} {
    }

    ForceField::~ForceField() = default;

    ForceField::ForceField(ForceField&& forceField)  noexcept {
        eps_inside_rc_ = forceField.eps_inside_rc_;
        eps_beyond_rc_ = forceField.eps_beyond_rc_;
        interactionsSpecs_ = std::move(forceField.interactionsSpecs_);
    }

    ForceField&
    ForceField::operator = (ForceField&& forceField) noexcept {
        eps_inside_rc_ = forceField.eps_inside_rc_;
        eps_beyond_rc_ = forceField.eps_beyond_rc_;
        interactionsSpecs_ = std::move(forceField.interactionsSpecs_);
        return *this;
    }

    void
    ForceField::addInteractionSpecification(const int_spec_t &spec) {
        interactionsSpecs_.emplace_back(spec);
    }

    const
    std::vector<ForceField::int_spec_t> &ForceField::interactionSpecifications() const {
        return interactionsSpecs_;
    }

    std::pair<real_t, real_t>
    ForceField::lennardJonesParameters(const spec_ptr_t &spec1,
                                       const spec_ptr_t &spec2) const
    {
        static util::Logger logger("simploce::ForceField::lennardJonesParameters()");

        for (auto& spec : this->interactionsSpecs_) {
            if ( spec.type == conf::LJ ) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.C12, spec.C6));
                }
            }
        }
        // Not found.
        logger.warn("(" + spec1->name() + ", " + spec2->name() +
                    "): No 'non-bonded Lennard-Jones' interaction parameters for this particle specifications pair.");
        return std::make_pair(0.0, 0.0);
    }

    std::pair<real_t, real_t>
    ForceField::harmonicParameters(const spec_ptr_t &spec1,
                                   const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::harmonicParameters()");
        for (auto& spec : this->interactionsSpecs_) {
            if ( spec.type == conf::HP ) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.r0, spec.fc));
                }
            }
        }
        // Not found.
        logger.warn("(" + spec1->name() + ", " + spec2->name() +
                    "): No 'harmonic bond' interaction parameters for this particle specifications pair.");
        return std::make_pair(0.0, 0.0);
    }

    std::pair<real_t, real_t>
    ForceField::halveAttractiveQuarticParameters(const spec_ptr_t &spec1,
                                                 const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::halveAttractiveQuarticParameters()");
        for (auto& spec : this->interactionsSpecs_) {
            if ( spec.type == conf::HA_QP ) {
                if ((spec.spec1 == spec1 && spec.spec2 == spec2) ||
                    (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.r0, spec.fc));
                }
            }
        }
        // Not found.
        logger.warn("(" + spec1->name() + ", " + spec2->name() +
                    "): No 'halve attractive quartic bond' interaction for this particle specifications pair.");
        return std::make_pair(0.0, 0.0);
    }

    std::tuple<real_t, real_t, real_t, real_t>
    ForceField::lennardJonesReactionFieldParameters(const spec_ptr_t &spec1,
                                                   const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::lennardJonesReactionFieldParameters()");
        for (auto& spec : this->interactionsSpecs_) {
            if ( spec.type == conf::LJ_RF ) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_tuple(spec.C12, spec.C6, eps_inside_rc_, eps_beyond_rc_));
                }
            }
        }
        // Not found.
        logger.warn("(" + spec1->name() + ", " + spec2->name() +
                    "): No 'non-bonded Lennard-Jones + Reaction Field' interaction parameters " +
                    "for this particle specifications pair.");
        return std::make_tuple(0.0, 0.0, 0.0, 0.0);

    }

    void
    ForceField::relativePermittivityInsideCutoff(real_t eps) {
        eps_inside_rc_ = eps;
    }

    real_t
    ForceField::relativePermittivityInsideCutoff() const {
        return eps_inside_rc_;
    }

    real_t
    ForceField::relativePermittivityBeyondCutoff() const {
        return eps_beyond_rc_;
    }

   void
   ForceField::relativePermittivityBeyondCutoff(real_t eps) {
        eps_beyond_rc_ = eps;
    }

    std::ostream&
    operator << (std::ostream& stream,
                 const ForceField& forceField) {
        stream.setf(std::ios::scientific);
        stream.precision(conf::PRECISION);
        const auto specs = forceField.interactionSpecifications();
        stream << "# Non-bonded" << std::endl;
        stream << "#     Type   Name #1   Name #2" << std::endl;
        for (auto& spec : specs) {
            if ( spec.type == conf::LJ ) {
                stream << std::setw(conf::NAME_WIDTH) << spec.type;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.C12;
                stream << std::setw(conf::REAL_WIDTH) << spec.C6;
                stream << conf::SPACE << "# C12 (kJ nm^12 mol^-1), C6 (kJ nm^6 mol^-1).";
                stream << std::endl;
            }
            // MORE HERE.
        }
        stream << "# Bonded" << std::endl;
        stream << "#     Type   Name #1   Name #2" << std::endl;
        for (auto& spec : specs) {
            if ( spec.type == conf::HP ) {
                stream << std::setw(conf::NAME_WIDTH) << spec.type;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.r0;
                stream << std::setw(conf::REAL_WIDTH) << spec.fc;
                stream << conf::SPACE << "# r0 (nm), fc (kJ/mol nm^2).";
                stream << std::endl;
            }
            // MORE HERE
        }
        return stream;
    }
}
