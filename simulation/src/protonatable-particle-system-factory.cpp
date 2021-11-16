/*
 * Author: André H. Juffer, Biocenter Oulu, University of Finland.
 *
 * Created on October 28, 2021.
 */

#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/util/logger.hpp"
#include <vector>

namespace simploce {

    ProtonatableParticleSystemFactory::ProtonatableParticleSystemFactory(const spec_catalog_ptr_t& catalog) :
            ParticleSystemFactory(catalog) {
    }

    prot_cg_sys
    ProtonatableParticleSystemFactory::protonatablePolarizableWater(const box_ptr_t& box,
                                                                    const density_t& densitySI,
                                                                    const temperature_t& temperature,
                                                                    std::size_t nLimit) {
        util::Logger logger("ProtonatableParticleSystemFactory::protonatablePolarizableWater");

        // Polarizable water.
        auto pWater = this->polarizableWater(box, densitySI, temperature, nLimit);

        // / Protonatable polarizable Water.
        prot_cg_sys ppWater;

        // Take beads in groups from pWater to create new beads in ppWater.
        auto catalog = this->catalog();
        pWater->doWithAllFreeGroups<void>([&ppWater, catalog, &logger] (
                std::vector<bead_ptr_t>& all,
                std::vector<bead_ptr_t>& free,
                std::vector<bead_group_ptr_t>& groups) {
            const std::string CW = "CW";
            const std::string DP = "DP";
            for (const auto& g : groups) {
                // Each water group consists of two beads, each of specifications CW and DP.
                auto beads = g->particles();
                std::vector<id_pair_t> bonds;
                std::vector<bead_ptr_t> newBeads;
                id_t idCW, idDP;
                for (const auto& bead : beads) {
                    auto spec = bead->spec();
                    if ( spec->name() == CW) {
                        auto sp = catalog->lookup("pCW");
                        std::string name = "pCW" + boost::lexical_cast<std::string>(all.size());
                        auto pCW = ppWater.addProtonatableBead(name, sp);
                        pCW->position(bead->position());
                        pCW->velocity(bead->velocity());
                        newBeads.emplace_back(pCW);
                        idCW = pCW->id();
                    } else if ( spec->name() == DP) {
                        auto sp = catalog->lookup("DP");
                        std::string name = "DP" + boost::lexical_cast<std::string>(all.size());
                        auto pDP = ppWater.addBead(name, sp);
                        pDP->position(bead->position());
                        pDP->velocity(bead->velocity());
                        newBeads.emplace_back(pDP);
                        idDP = pDP->id();
                    } else {
                        util::logAndThrow(logger, "Required specifications for CW and/or DP missing.");
                    }
                }
                bonds.emplace_back(std::make_pair(idCW, idDP));
                ppWater.addBeadGroup(newBeads, bonds);
            }
        });
        ppWater.box(box);

        logger.info("Number of protonatable beads: " +
                    boost::lexical_cast<std::string>(ppWater.numberProtonatableBeads()));

        return std::move(ppWater);
    }

    prot_cg_sys
    ProtonatableParticleSystemFactory::formicAcid(const box_ptr_t& box,
                                                  const density_t& densitySI,
                                                  const molarity_t& molarity,
                                                  const temperature_t& temperature) {

        util::Logger logger{"ParticleModelFactory::formicAcid"};
        const auto catalog = this->catalog();

        logger.info("Creating coarse grained particle model for formic acid (HCOOH) in water");

        logger.info("Requested molarity: " + boost::lexical_cast<std::string>(molarity));

        // Create CG protonatable polarizable water particle model.
        auto ppModel = this->protonatablePolarizableWater(box, densitySI, temperature);

        // Molarity to number of HCOOH molecules.
        volume_t volume = box->volume();
        number_density_t numberDensity =
            molarity() * units::si<real_t>::NA / units::mu<real_t>::l_to_nm3;
        int nHCOOH = util::nint(numberDensity() * volume());
        logger.info("Requested number of HCOOH: " + boost::lexical_cast<std::string>(nHCOOH));
        logger.debug("Replacing " + boost::lexical_cast<std::string>(nHCOOH) + " water groups by HCOOH beads.");
        auto spec = this->catalog()->lookup("HCOOH");
        logger.debug("HCOOH specification:");
        util::log(logger, *spec);
        ppModel.replaceGroupsByParticles(spec, nHCOOH);
        logger.info("Number of HCOOH: " + boost::lexical_cast<std::string>(ppModel.numberOfSpecifications(spec)));
        logger.info("Number of pCW: " + boost::lexical_cast<std::string>(ppModel.numberOfSpecifications(catalog->lookup("pCW"))));
        return std::move(ppModel);
    }
}
