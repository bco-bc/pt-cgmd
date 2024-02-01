/*
 * Created by André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#include "simploce/potentials/force-field.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/logger.hpp"
#include <boost/algorithm/string.hpp>
#include <utility>
#include <vector>
#include <memory>

namespace simploce {

    static void notFound(util::Logger& logger,
                         const spec_ptr_t& spec1, const spec_ptr_t& spec2) {
        util::logAndThrow(logger,
                          "(" + spec1->name() + ", " + spec2->name() +
                          "): No interaction parameters for this particle specifications pair.");
    }

    ff_ptr_t
    ForceField::obtainFrom(std::istream &stream,
                           const spec_catalog_ptr_t &catalog) {
        util::Logger logger("simploce::ForceField::obtainFrom()");
        logger.trace("Entering.");

        ff_ptr_t forceField = std::make_shared<ForceField>();
        char charBuffer[100];
        std::string stringBuffer;

        // Skip two headers.
        std::getline(stream, stringBuffer);
        std::getline(stream, stringBuffer);

        // Read non-bonded interaction parameters.
        int nNonBonded;
        stream >> nNonBonded;
        logger.debug(std::to_string(nNonBonded) + " : Number of non-bonded interaction types.");
        std::getline(stream, stringBuffer);
        for (int k = 0; k != nNonBonded; ++k) {
            // Read interaction specification type name.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string typeName = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(typeName);
            logger.debug(typeName + ": Got non-bonded pair potential.");

            // Get two particle specifications.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName1 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName1);
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string specName2 = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(specName2);
            std::string msg = specName1;
            msg += ", ";
            msg += specName2 + ": Particle specification pair.";
            logger.debug(msg);
            auto spec1 = catalog->lookup(specName1);
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
            } else if (typeName == conf::SR) {
                spec.typeName = conf::SR;
                stream >> spec.max_a;
            } else if (typeName == conf::GAUSS_SF) {
                spec.typeName = conf::GAUSS_SF;
                stream >> spec.sigma1 >> spec.sigma2;
            } else if (typeName == conf::GAUSS_SF_SR) {
                spec.typeName = conf::GAUSS_SF_SR;
                stream >> spec.sigma1 >> spec.sigma2 >> spec.max_a;
            } else if (typeName == conf::NONE_INTERACTING) {
                spec.typeName = conf::NONE_INTERACTING;
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
        logger.debug(std::to_string(nBonded) + ": Number of bonded interaction types.");
        std::getline(stream, stringBuffer);
        for (int k = 0; k != nBonded; ++k) {
            // Read interaction specification type name.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string typeName = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(typeName);
            logger.debug(typeName + ": Got bonded pair potential.");

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
            if (typeName == conf::HP || typeName == conf::HA_QP || typeName == conf::HA_HP) {
                spec.typeName = typeName;
                stream >> spec.r0 >> spec.fc;
            } else if (typeName == conf::SR) {
                spec.typeName = conf::SR;
                stream >> spec.max_a;
            } else if (typeName == conf::HP_SR) {
                spec.typeName = typeName;
                stream >> spec.r0 >> spec.fc >> spec.max_a;
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
        logger.debug(std::to_string(nExternal) + ": Number of external potential types.");
        std::getline(stream, stringBuffer);
        for (int k = 0; k != nExternal; ++k) {
            // Read interaction specification type name.
            stream.read(charBuffer, conf::NAME_WIDTH);
            std::string typeName = std::string(charBuffer, conf::NAME_WIDTH);
            boost::trim(typeName);
            logger.debug(typeName + ": Got external potential.");

            int_spec_t spec;
            if (typeName == conf::VOLTAGE) {
                spec.typeName = conf::VOLTAGE;
                real_t ex, ey, ez;
                stream >> ex >> ey >> ez >> spec.eps_r;
                spec.e0 = el_field_t{ex, ey, ez};
            } else if (typeName == conf::WALL) {
                spec.typeName = conf::WALL;
                std::string s;
                stream >> spec.C12 >> spec.C6 >> spec.sigma >> spec.distance >> s;
                boost::trim(s);
                spec.plane = s;
            } else if (typeName == conf::PRESSURE_GRADIENT) {
                spec.typeName = conf::PRESSURE_GRADIENT;
                real_t fx, fy, fz;
                stream >> fx >> fy >> fz;
                spec.fe = force_t(fx, fy, fz);
            } else {
                util::logAndThrow(logger, "'" + typeName + "': No such external potential.");
            }

            forceField->addExternalSpecification(spec);
            std::getline(stream, stringBuffer);
        }

        // Log some information.
        std::size_t numberNonBonded = forceField->nonBondedSpecifications().size();
        std::size_t numberBonded = forceField->bondedSpecifications().size();
        logger.info(std::to_string(numberNonBonded) + ": Number of non-bonded interaction potential types.");
        logger.info(std::to_string(numberBonded) + ": Number of bonded interaction potential types.");
        logger.info(std::to_string(nExternal) + ": Number of external potential potential types.");

        logger.trace("Leaving.");
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
        notFound(logger, spec1, spec2);
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
        notFound(logger, spec1, spec2);
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
        notFound(logger, spec1, spec2);
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
        notFound(logger, spec1, spec2);
        return 0.0;
    }

    real_t
    ForceField::hardSphereShiftedForce(const spec_ptr_t &spec1,
                                       const spec_ptr_t &spec2) const {
        return this->shiftedForce(spec1, spec2);
    }

    std::pair<real_t, real_t>
    ForceField::gaussianChargeDensity(const simploce::spec_ptr_t &spec1, const simploce::spec_ptr_t &spec2) const {
        util::Logger logger("simploce::ForceField::gaussianChargeDensity()");
        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::GAUSS_SF) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.sigma1, spec.sigma2));
                }
            }
        }
        // Not found.
        notFound(logger, spec1, spec2);
        return std::make_pair(0.0, 0.0);
    }

    std::tuple<real_t, real_t, real_t>
    ForceField::gaussianChargeDensitySoftRepulsion(const simploce::spec_ptr_t &spec1,
                                                   const simploce::spec_ptr_t &spec2) const {
        util::Logger logger("simploce::ForceField::gaussianChargeDensitySoftRepulsion()");
        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::GAUSS_SF_SR) {
                if ( (spec.spec1 == spec1 && spec.spec2 == spec2) ||
                     (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_tuple(spec.sigma1, spec.sigma2, spec.max_a));
                }
            }
        }
        // Not found
        notFound(logger, spec1, spec2);
        return std::make_tuple(0.0, 0.0, 0.0);
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
        notFound(logger, spec1, spec2);
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
        notFound(logger, spec1, spec2);
        return std::make_pair(0.0, 0.0);
    }

    std::pair<real_t, real_t>
    ForceField::halveAttractiveHarmonic(const spec_ptr_t &spec1, const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::halveAttractiveHarmonic()");
        for (auto& spec : this->bondedSpecs_) {
            if (spec.typeName == conf::HA_HP ) {
                if ((spec.spec1 == spec1 && spec.spec2 == spec2) ||
                    (spec.spec1 == spec2 && spec.spec2 == spec1) ) {
                    return std::move(std::make_pair(spec.r0, spec.fc));
                }
            }
        }
        // Not found.
        notFound(logger, spec1, spec2);
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
        notFound(logger, spec1, spec2);
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
        notFound(logger, spec1, spec2);
        return std::move(std::make_tuple(0.0, 0.0, 0.0));
    }

    real_t
    ForceField::softRepulsion(const spec_ptr_t &spec1,
                              const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::softRepulsion()");
        for (auto& spec : this->nonBondedSpecs_) {
            if (spec.typeName == conf::SR ) {
                if ( (*(spec.spec1) == *spec1 && *(spec.spec2) == *spec2) ||
                     (*(spec.spec1) == *spec2 && *(spec.spec2) == *spec1) ) {
                    return spec.max_a;
                }
            }
        }
        // Not found.
        notFound(logger, spec1, spec2);
        return 0.0;
    }

    std::tuple<real_t, real_t, real_t>
    ForceField::harmonicSoftRepulsion(const spec_ptr_t &spec1,
                                      const spec_ptr_t &spec2) const {
        static util::Logger logger("simploce::ForceField::harmonicSoftRepulsion()");
        for (auto& spec : this->bondedSpecs_) {
            if (spec.typeName == conf::HP_SR ) {
                if ( (*(spec.spec1) == *spec1 && *(spec.spec2) == *spec2) ||
                     (*(spec.spec1) == *spec2 && *(spec.spec2) == *spec1) ) {
                    return std::make_tuple(spec.r0, spec.fc, spec.max_a);
                }
            }
        }

        // Not found.
        notFound(logger, spec1, spec2);
        return std::make_tuple(0.0, 0.0, 0.0);
    }

    std::pair<srf_charge_density_t, real_t>
    ForceField::constSurfaceChargeDensity() const {
        return std::pair<srf_charge_density_t, real_t>{0.0, 0.0};
    }

    std::tuple<el_pot_diff_t, dist_t, real_t>
    ForceField::electricPotentialDifference() const {
        return std::tuple<el_pot_diff_t, dist_t, real_t>{0.0, 0.0, 0.0};
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
            } else if (spec.typeName == conf::SR) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.max_a;
                stream << conf::SPACE << "# Maximum repulsion.";
                stream << std::endl;
            } else if (spec.typeName == conf::GAUSS_SF) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.sigma1;
                stream << std::setw(conf::REAL_WIDTH) << spec.sigma2;
                stream << conf::SPACE << "# Widths of charge densities";
                stream << std::endl;
            } else if (spec.typeName == conf::GAUSS_SF_SR) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.sigma1;
                stream << std::setw(conf::REAL_WIDTH) << spec.sigma2;
                stream << std::setw(conf::REAL_WIDTH) << spec.max_a;
                stream << conf::SPACE << "# Widths of charge densities, maximum repulsion).";
                stream << std::endl;
            } else if (spec.typeName == conf::HS_SF) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.eps_inside_rc;
                stream << conf::SPACE << "# eps_inside_rc";
            } else if (spec.typeName == conf::NONE_INTERACTING) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << conf::SPACE << "# No interaction.";
                stream << std::endl;
            }
        }
        const auto& bSpecs = forceField.bondedSpecifications();
        stream << "# Bonded" << std::endl;
        stream << "#     Type   Name #1   Name #2" << std::endl;
        stream << bSpecs.size() << std::endl;
        for (auto& spec : bSpecs) {
            if (spec.typeName == conf::HP || spec.typeName == conf::HA_QP || spec.typeName == conf::HA_HP ) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.r0;
                stream << std::setw(conf::REAL_WIDTH) << spec.fc;
                stream << conf::SPACE << "# r0, fc";
                stream << std::endl;
            } else if (spec.typeName == conf::SR) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.max_a;
                stream << conf::SPACE << "# Maximum repulsion (DPD).";
                stream << std::endl;
            } else if (spec.typeName == conf::HP_SR) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::NAME_WIDTH) << spec.spec1->name();
                stream << std::setw(conf::NAME_WIDTH) << spec.spec2->name();
                stream << std::setw(conf::REAL_WIDTH) << spec.r0;
                stream << std::setw(conf::REAL_WIDTH) << spec.fc;
                stream << std::setw(conf::REAL_WIDTH) << spec.max_a;
                stream << conf::SPACE << "# r0, fc, aij (maximum repulsion, DPD)";
                stream << std::endl;
            }
        }
        const auto& exSpecs = forceField.externalSpecifications();
        stream << "# External" << std::endl;
        stream << "#     Type    Parameters" << std::endl;
        stream << exSpecs.size() << std::endl;
        for (auto& spec: exSpecs) {
            if (spec.typeName == conf::VOLTAGE) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << std::setw(conf::REAL_WIDTH) << spec.e0;
                stream << conf::SPACE << "# z-component of external static electric homogeneous electric field.";
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
            } else if (spec.typeName == conf::PRESSURE_GRADIENT) {
                stream << std::setw(conf::NAME_WIDTH) << spec.typeName;
                stream << spec.fe << "# Component of external force  (constant, homogeneous).";
                stream << std::endl;
            }
        }
        return stream;
    }
}
