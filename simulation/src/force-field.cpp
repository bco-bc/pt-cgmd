/*
 * Created by André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#include "simploce/potentials/force-field.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/types/s-types.hpp"
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
        util::Logger logger("simploce::ForceField::forceField");

        ff_ptr_t forceField = std::make_shared<ForceField>();
        char charBuffer[100];
        std::string stringBuffer;

        // Skip two headers.
        std::getline(stream, stringBuffer);
        std::getline(stream, stringBuffer);

        // Read non-bonded interaction parameters.
        int nNonBonded;
        stream >> nNonBonded;
        std::getline(stream, stringBuffer);
        for (int k = 0; k != nNonBonded; ++k) {
            // Read interaction specification type name.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string typeName = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(typeName);

            // Get two particle specifications.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName1 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName1);
            auto spec1 = catalog->lookup(specName1);
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName2 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName2);
            auto spec2 = catalog->lookup(specName2);

            // Read interaction parameters.
            int_spec_t spec;
            spec.spec1 = spec1;
            spec.spec2 = spec2;
            if (typeName == conf::LJ) {
                spec.typeName = typeName;
                stream >> spec.C12 >> spec.C6;
            } else if (typeName == conf::LJ_RF) {
                spec.typeName = typeName;
                stream >> spec.C12 >> spec.C6 >> spec.eps_inside_rc >> spec.eps_outside_rc;
            } else if (typeName == conf::RF) {
                spec.typeName = typeName;
                stream >> spec.eps_inside_rc >> spec.eps_outside_rc;
            } else if (typeName == conf::HS_SF || typeName == conf::SF ||
                       typeName == conf::HS_SC || typeName == conf::SC) {
                spec.typeName = typeName;
                stream >> spec.eps_inside_rc;
            } else if (typeName == conf::LJ_SF) {
                spec.typeName = conf::LJ_SF;
                stream >> spec.C12 >> spec.C6 >> spec.eps_inside_rc;
            } else {
                util::logAndThrow(logger, "'" + typeName + "': No such interaction type.");
            }
            forceField->addNonBondedSpecification(spec);

            // Read and ignore description.
            std::getline(stream, stringBuffer);
        }

        // Skip two headers.
        std::getline(stream, stringBuffer);
        std::getline(stream, stringBuffer);

        int nBonded;
        stream >> nBonded;
        std::getline(stream, stringBuffer);
        for (int k = 0; k != nBonded; ++k) {
            // Read interaction specification type name.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string typeName = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(typeName);

            // Get two particle specifications.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName1 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName1);
            auto spec1 = catalog->lookup(specName1);
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName2 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName2);
            auto spec2 = catalog->lookup(specName2);

            int_spec_t spec;
            spec.spec1 = spec1;
            spec.spec2 = spec2;
            if (typeName == conf::HP || typeName == conf::HA_QP ) {
                spec.typeName = typeName;
                stream >> spec.r0 >> spec.fc;
            } else {
                util::logAndThrow(logger, "'" + typeName + "': No such interaction type.");
            }
            forceField->addBondedSpecification(spec);

            // Read and ignore description.
            std::getline(stream, stringBuffer);
        }

        // Skip two headers.
        std::getline(stream, stringBuffer);
        std::getline(stream, stringBuffer);

        int nExternal;
        stream >> nExternal;
        std::getline(stream, stringBuffer);
        for (int k = 0; k != nExternal; ++k) {
            // Read interaction specification type name.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string typeName = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(typeName);

            int_spec_t spec;
            if (typeName == conf::ELECTRIC_POTENTIAL_DIFFERENCE) {
                spec.typeName = conf::ELECTRIC_POTENTIAL_DIFFERENCE;
                std::string s;
                stream >> spec.deltaV >> spec.distance >> spec.eps_inside_rc >> s;
                boost::trim(s);
                spec.direction = s[0];
            } else if (typeName == conf::WALL) {
                spec.typeName = conf::WALL;
                std::string s;
                stream >> spec.C12 >> spec.C6 >> spec.sigma >> spec.distance >> s;
                boost::trim(s);
                spec.plane = s;
            } else {
                util::logAndThrow(logger, "'" + typeName + "': No such external potential.");
            }

            forceField->addExternalSpecification(spec);
            std::getline(stream, stringBuffer);
        }

        // Log some information.
        std::size_t numberNonBonded = forceField->nonBondedSpecifications().size();
        std::size_t numberBonded = forceField->bondedSpecifications().size();
        logger.debug("Number of non-bonded interaction potential types: " +
                     util::toString(numberNonBonded));
        logger.debug("Number of bonded interaction potential types: " +
                     util::toString(numberBonded));
        logger.debug("Number of external potential potential types: " +
                     util::toString(nExternal));

        return std::move(forceField);
    }

    ForceField::ForceField() : bondedSpecs_{}, nonBondedSpecs_{} {
    }

    ForceField::~ForceField() = default;

    ForceField::ForceField(ForceField&& forceField)  noexcept {
        bondedSpecs_ = std::move(forceField.bondedSpecs_);
        nonBondedSpecs_ = std::move(forceField.nonBondedSpecs_);
    }

    ForceField&
    ForceField::operator = (ForceField&& forceField) noexcept {
        bondedSpecs_ = std::move(forceField.bondedSpecs_);
        nonBondedSpecs_ = std::move(forceField.nonBondedSpecs_);
        return *this;
    }

    void
    ForceField::addNonBondedSpecification(const int_spec_t &spec) {
        nonBondedSpecs_.emplace_back(spec);
    }

    void
    ForceField::addBondedSpecification(const int_spec_t &spec) {
        bondedSpecs_.emplace_back(spec);
    }

    void
    ForceField::addExternalSpecification(const int_spec_t& spec) {
        external_.emplace_back(spec);
    }

    const std::vector<ForceField::int_spec_t>&
    ForceField::nonBondedSpecifications() const {
        return nonBondedSpecs_;
    }

    const std::vector<ForceField::int_spec_t>&
    ForceField::bondedSpecifications() const {
        return bondedSpecs_;
    }

    const std::vector<ForceField::int_spec_t>&
    ForceField::externalSpecifications() const {
        return external_;
    }

    std::pair<real_t, real_t>
    ForceField::lennardJones(const spec_ptr_t &spec1,
                             const spec_ptr_t &spec2) const
    {
        static util::Logger logger("simploce::ForceField::lennardJones()");

        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::LJ ) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.C12, spec.C6));
                }
            }
        }
        // Not found.
        util::logAndThrow(logger,
                          "(" + spec1->name() + ", " + spec2->name() +
                          "): No interaction parameters for this particle specifications pair.");
        return std::move(std::make_pair(0.0, 0.0));
    }

    std::pair<real_t, real_t>
    ForceField::reactionField(const spec_ptr_t &spec1,
                              const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::reactionField()");

        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::RF || spec.typeName == conf::HS_RF) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.eps_inside_rc, spec.eps_outside_rc));
                }
            }
        }
        // Not found.
        util::logAndThrow(logger,
                          "(" + spec1->name() + ", " + spec2->name() +
                          "): No interaction parameters for this particle specifications pair.");
        return std::move(std::make_pair(0.0, 0.0));
    }

    std::pair<real_t, real_t>
    ForceField::hardSphereReactionField(const spec_ptr_t &spec1,
                                        const spec_ptr_t &spec2) const {
        return std::move(this->reactionField(spec1, spec2));
    }

    real_t
    ForceField::screenedCoulomb(const spec_ptr_t &spec1,
                                const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::screenedCoulomb()");
        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::SC || spec.typeName == conf::HS_SC) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return spec.eps_inside_rc;
                }
            }
        }
        // Not found.
        util::logAndThrow(logger,
                          "(" + spec1->name() + ", " + spec2->name() +
                          "): No interaction parameters for this particle specifications pair.");
        return 0.0;
    }

    real_t
    ForceField::hardSphereScreenedCoulomb(const spec_ptr_t &spec1,
                                          const spec_ptr_t &spec2) const {
        return this->screenedCoulomb(spec1, spec2);
    }

    real_t
    ForceField::shiftedForce(const spec_ptr_t &spec1,
                             const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::shiftedForce()");
        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::SF || spec.typeName == conf::HS_SF) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return spec.eps_inside_rc;
                }
            }
        }
        // Not found.
        util::logAndThrow(logger,
                          "(" + spec1->name() + ", " + spec2->name() +
                          "): No interaction parameters for this particle specifications pair.");
        return 0.0;
    }

    real_t
    ForceField::hardSphereShiftedForce(const spec_ptr_t &spec1,
                                       const spec_ptr_t &spec2) const {
        return this->shiftedForce(spec1, spec2);
    }

    std::pair<real_t, real_t>
    ForceField::harmonic(const spec_ptr_t &spec1,
                         const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::harmonic()");
        for (auto& spec : this->bondedSpecs_) {
            if (spec.typeName == conf::HP ) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.r0, spec.fc));
                }
            }
        }
        // Not found.
        util::logAndThrow(logger,
                          "(" + spec1->name() + ", " + spec2->name() +
                          "): No interaction parameters for this particle specifications pair.");
        return std::make_pair(0.0, 0.0);
    }

    std::pair<real_t, real_t>
    ForceField::halveAttractiveQuartic(const spec_ptr_t &spec1,
                                       const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::halveAttractiveQuartic()");
        for (auto& spec : this->bondedSpecs_) {
            if (spec.typeName == conf::HA_QP ) {
                if ((spec.spec1 == spec1 && spec.spec2 == spec2) ||
                    (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.r0, spec.fc));
                }
            }
        }
        // Not found.
        logger.warn("(" + spec1->name() + ", " + spec2->name() +
                    "): No interaction for this particle specifications pair.");
        return std::make_pair(0.0, 0.0);
    }

    std::tuple<real_t, real_t, real_t, real_t>
    ForceField::lennardJonesReactionField(const spec_ptr_t &spec1,
                                          const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::lennardJonesReactionField()");
        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::LJ_RF ) {
                if ( (*(spec.spec1) == *spec1 && *(spec.spec2) == *spec2) ||
                     (*(spec.spec1) == *spec2 && *(spec.spec2) == *spec1) ) {
                    return std::move(std::make_tuple(spec.C12, spec.C6, spec.eps_inside_rc, spec.eps_outside_rc));
                }
            }
        }
        // Not found.
        logger.warn("(" + spec1->name() + ", " + spec2->name() +
                    "): No interaction parameters for this particle specifications pair.");
        return std::move(std::make_tuple(0.0, 0.0, 0.0, 0.0));
    }

    std::tuple<real_t, real_t, real_t>
    ForceField::lennardJonesShiftedForce(const spec_ptr_t &spec1,
                                         const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::lennardJonesShiftedForce()");
        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::LJ_SF ) {
                if ( (*(spec.spec1) == *spec1 && *(spec.spec2) == *spec2) ||
                     (*(spec.spec1) == *spec2 && *(spec.spec2) == *spec1) ) {
                    return std::move(std::make_tuple(spec.C12, spec.C6, spec.eps_inside_rc));
                }
            }
        }
        // Not found.
        logger.warn("(" + spec1->name() + ", " + spec2->name() +
                    "): No interaction parameters for this particle specifications pair.");
        return std::move(std::make_tuple(0.0, 0.0, 0.0));
    }

    std::pair<srf_charge_density_t, real_t>
    ForceField::constSurfaceChargeDensity() const {
        return std::pair<srf_charge_density_t, real_t>{0.0, 0.0};
    }

    std::tuple<el_pot_diff, dist_t, real_t>
    ForceField::electricPotentialDifference() const {
        return std::tuple<el_pot_diff, dist_t, real_t>{0.0, 0.0, 0.0};
    }

    std::ostream&
    operator << (std::ostream& stream, const ForceField& forceField) {
        stream.setf(std::ios::scientific);
        stream.precision(conf::PRECISION);
        const auto& nbSpecs = forceField.nonBondedSpecifications();
        stream << "# Non-bonded" << std::endl;
        stream << "#     Type   Name #1   Name #2" << std::endl;
        stream << nbSpecs.size() << std::endl;
        for (auto& spec : nbSpecs) {
            if (spec.typeName == conf::LJ ) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.C12;
                stream << std::setw(conf::REAL_WIDTH) << spec.C6;
                stream << conf::SPACE << "# C12, C6";
                stream << std::endl;
            } else if (spec.typeName == conf::LJ_RF) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.C12;
                stream << std::setw(conf::REAL_WIDTH) << spec.C6;
                stream << std::setw(conf::REAL_WIDTH) << spec.eps_inside_rc;
                stream << std::setw(conf::REAL_WIDTH) << spec.eps_outside_rc;
                stream << conf::SPACE << "# C12, C6, eps_inside_rc, eps_outside_rc";
                stream << std::endl;
            } else if (spec.typeName == conf::LJ_SF) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.C12;
                stream << std::setw(conf::REAL_WIDTH) << spec.C6;
                stream << std::setw(conf::REAL_WIDTH) << spec.eps_inside_rc;
                stream << conf::SPACE << "# C12, C6, eps_inside_rc";
                stream << std::endl;
            } else if (spec.typeName == conf::SC || spec.typeName == conf::HS_SC) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.eps_inside_rc;
                stream << conf::SPACE << "# eps_inside_rc";
                stream << std::endl;
            } else if (spec.typeName == conf::RF || spec.typeName == conf::HS_RF) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.eps_inside_rc;
                stream << std::setw(conf::REAL_WIDTH) << spec.eps_outside_rc;
                stream << conf::SPACE << "# eps_inside_rc, eps_outside_rc";
                stream << std::endl;
            }
        }
        const auto& bSpecs = forceField.bondedSpecifications();
        stream << "# Bonded" << std::endl;
        stream << "#     Type   Name #1   Name #2" << std::endl;
        stream << bSpecs.size() << std::endl;
        for (auto& spec : bSpecs) {
            if (spec.typeName == conf::HP || spec.typeName == conf::HA_QP ) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.r0;
                stream << std::setw(conf::REAL_WIDTH) << spec.fc;
                stream << conf::SPACE << "# r0, fc";
                stream << std::endl;
            }
        }
        const auto& exSpecs = forceField.externalSpecifications();
        stream << "# External" << std::endl;
        stream << "#     Type    Parameters" << std::endl;
        stream << exSpecs.size() << std::endl;
        for (auto& spec: exSpecs) {
            if (spec.typeName == conf::ELECTRIC_POTENTIAL_DIFFERENCE) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::REAL_WIDTH) << spec.deltaV;
                stream << std::setw(conf::REAL_WIDTH) << spec.distance;
                stream << conf::SPACE << spec.eps_inside_rc;
                stream << conf::SPACE << spec.direction;
                stream << conf::SPACE << "# deltaV distance eps_r direction";
                stream << std::endl;
            } else if (spec.typeName == conf::WALL) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << conf::SPACE << spec.C12;
                stream << conf::SPACE << spec.C6;
                stream << conf::SPACE << spec.sigma;
                stream << conf::SPACE << spec.distance;
                stream << conf::SPACE << spec.plane;
                stream << conf::SPACE << "# C12, C6, sigma, distance, plane";
                stream << std::endl;
            }
        }
        return stream;
    }
}
