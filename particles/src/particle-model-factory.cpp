/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   particle-model-factory.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 3, 2019, 12:50 PM
 */

#include "simploce/particle/particle-model-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/polarizable-water.hpp"
#include "simploce/particle/pfactory.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/continuous-protonatable-bead.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/mu-units.hpp"
#include "simploce/util/box.hpp"
#include <stdexcept>
#include <iostream>

namespace simploce {
    
    static const length_t DISTANCE_CW_DP = 0.2; // nm.
    
    // From Maxwell velocity distribution.
    // https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    template <typename P>
    void assignMomentum(std::shared_ptr<P>& particle, const temperature_t& temperature)
    {
        static real_t KB = MUUnits<real_t>::KB;
        static std::random_device rd{};
        static std::mt19937 gen{rd()};
        static bool init = false;
        if ( !init ) {
            gen.seed(util::seedValue<std::size_t>());  // seed() From util.hpp
            init = true;
        }

        // For each component, take its value from a Gaussian density for velocity with
        // standard deviation 'sigma' and zero average.
        real_t sigma = std::sqrt(KB * temperature / particle->mass());
        std::normal_distribution<real_t> gaussian{0.0, sigma};  // Zero average.
        momentum_t p{};
            for (std::size_t k = 0; k != 3; ++k) {
                real_t v = gaussian(gen);
            p[k] = particle->mass()() * v;
        }
        particle->momentum(p);
    }
    
    
    ParticleModelFactory::ParticleModelFactory(const spec_catalog_ptr_t& catalog) :
        catalog_{catalog}
    {        
        if ( !catalog_ ) {
            throw std::domain_error("ParticleModelFactory: Missing particle specifications catalog.");
        }
    }
        
    cg_pol_water_ptr_t
    ParticleModelFactory::polarizableWater(const box_ptr_t& box,
                                           const density_t atDensitySI,
                                           const temperature_t temperature,
                                           bool protonatable) const
    {
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        std::clog << "Creating coarse grained particle model "
                     "for polarizable water." << std::endl;
        
        std::clog << "Temperature: " << temperature << std::endl;

        // Convert kg/m^3 to u/nm^3.
        density_t atDensity = atDensitySI / (SIUnits<real_t>::MU * 1.0e+27);
        std::clog << "Requested density (kg/m^3): " << atDensitySI << std::endl;
        std::clog << "Requested density (u/nm^3): " << atDensity << std::endl;
        
        // Spacing between particles.
        const length_t spacing{0.53};  // Roughly the location of the first peak of g(r), in nm.
        std::clog << "Spacing between particles: " << spacing << " nm" << std::endl;
        
        // Box details.
        length_t Lx = box->lengthX();
        length_t Ly = box->lengthY();
        length_t Lz = box->lengthZ();
        volume_t volume = box->volume();
        std::clog << "Box size (nm): " << box->size() << std::endl;
        std::clog << "Box volume (nm^3): " << volume << std::endl;
        
        // How many CG particles? A single polarizable CG water particle represents 
        // 5 water molecules. Each CG water particles consists of two connected 
        // CG particles (CW and DP).
        spec_ptr_t mh2o = catalog_->molecularWater();
        spec_ptr_t cwSpec = 
            protonatable ? catalog_->lookup("PCW") : catalog_->lookup("CW");
        spec_ptr_t dpSpec = catalog_->lookup("DP");
        std::clog << "\"Ideal\" distance between CW and DP: " << DISTANCE_CW_DP
                  << " nm" << std::endl;

        number_density_t atNumberDensity = atDensity / mh2o->mass();
        std::size_t natWaters = atNumberDensity * volume;
        //std::size_t ncgWaters = util::nint(real_t(natWaters) / 5.0);
        std::size_t ncgWaters = real_t(natWaters) / 5.0;
        number_density_t cgNumberDensity = real_t(ncgWaters) / volume();
        density_t cgDensity = cgNumberDensity * (cwSpec->mass() + dpSpec->mass());
        std::clog << "CG: Requested density (u/nm^3): " << cgDensity << std::endl;
        std::clog << "AT: Requested number density (1/nm^3): " << atNumberDensity << std::endl;
        std::clog << "AT: Requested number density (1/m^3): " << atNumberDensity*1.0e+27 << std::endl;
        std::clog << "CG: Requested number density (1/nm^3): " << cgNumberDensity << std::endl;
        std::clog << "CG: Requested number density (1/m^3): " << cgNumberDensity*1.0e+27 << std::endl;
        std::clog << "AT: Requested number of water molecules: " << natWaters << std::endl;
        std::clog << "CG: Requested number of water particles: " << ncgWaters << std::endl;
        std::clog.flush();

        std::size_t nx = util::nint(Lx / spacing);
        std::size_t ny = util::nint(Ly / spacing);
        std::size_t nz = util::nint(Lz / spacing);
        std::clog << "Distance spacing between CG water particles: " << spacing << std::endl;
        std::clog << "Number of coordinates in x-direction: " << nx << std::endl;
        std::clog << "Number of coordinates in y-direction: " << ny << std::endl;
        std::clog << "Number of coordinates in z-direction: " << nz << std::endl;
        std::clog.flush();
        
        // Start from an empty particle model.
        cg_pol_water_ptr_t cg = std::make_shared<PolarizableWater>();
        
        // Add beads.
        std::size_t counter = 0;
        std::size_t i = 0, j = 0, k = 0;         
        while ( i < nx && counter < ncgWaters ) {
            while ( j < ny && counter < ncgWaters ) {
                while ( k < nz && counter < ncgWaters ) {
                    
                    std::size_t id = 2 * counter + 1;

                    // Single water group.
                    std::vector<bead_ptr_t> beads{};
                    std::vector<id_pair_t> bbonds{};
                            
                    real_t x = (i + 0.5) * spacing();
                    real_t y = (j + 0.5) * spacing();
                    real_t z = (k + 0.5) * spacing();
                    position_t r1{x,y,z};
                    bead_ptr_t cwBead = 
                        protonatable ? 
                        cg->addContinuousProtonatableBead(id, "P-CW", r1, 0, cwSpec, false) :
                        cg->addBead(id, "CW", r1, cwSpec, false);
                    assignMomentum(cwBead, temperature);
                    beads.push_back(cwBead);
          
                    // Place the DP parallel to one of the 3 coordinate axes.
                    position_t r2{};
                    std::size_t l = util::random<real_t>() * 3.0;
                    switch (l) {
                        case 0: {
                            real_t coord = x + DISTANCE_CW_DP();
                            r2 = position_t{coord, y, z};              
                            break;
                        }
                        case 1: {
                            real_t coord = y + DISTANCE_CW_DP();
                            r2 = position_t{x,coord,z};
                            break;
                        }
                        default: {
                            real_t coord = z + DISTANCE_CW_DP();
                            r2 = position_t{x,y,coord};
                            break;
                        }
                    }
                    id += 1;
                    bead_ptr_t dpBead = cg->addBead(id, "DP", r2, dpSpec, false);
                    assignMomentum<Bead>(dpBead, temperature);
                    beads.push_back(dpBead);
                    
                    auto pair = std::make_pair<int, int>(cwBead->id(), dpBead->id());
                    bbonds.push_back(pair);
                    
                    cg->addBeadGroup(beads, bbonds);
          
                    counter += 1;
                    k += 1;
                }
                k = 0;
                j += 1;
            }
            j = 0;
            i += 1;
        }
        std::clog << "Created " << cg->numberOfParticles() << " beads." << std::endl;
        std::clog << "Created " << cg->numberOfParticleGroups() << " waters groups." << std::endl;
        
        return cg;
    }
    
    cg_ptr_t 
    ParticleModelFactory::formicAcidSolution(const box_ptr_t& box,
                                             const density_t atDensitySI,
                                             const molarity_t molarity,
                                             const temperature_t temperature,
                                             bool protonatable) const
    {
        std::clog << "Creating coarse grained particle model for formic acid "
                     "(HCOOH) in water." << std::endl;

        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        std::clog << "Molarity: " << molarity << " M" << std::endl;
        
        // Create CG polarizable water.
        cg_pol_water_ptr_t cg = 
            this->polarizableWater(box, atDensitySI, temperature, protonatable);
        
        // Molarity to number of HCOOH molecules.
        volume_t volume = box->volume();
        number_density_t numberDensity =
            molarity() * SIUnits<real_t>::NA / MUUnits<real_t>::l_to_nm3;
        std::size_t nHCOOH = util::nint<real_t>(numberDensity() * volume());
        std::clog << "Number of HCOOH: " << nHCOOH << std::endl;
        
        std::clog << "Replacing " << nHCOOH << " water groups by HCOOH beads." << std::endl;
        spec_ptr_t spec = catalog_->lookup("HCOOH");
        auto size = cg->numberOfBeads();
        for (std::size_t counter = 0; counter != nHCOOH; ++counter) {
            std::size_t id = size + counter + 1;
            position_t r = cg->removeGroup_();
            auto bead = cg->addContinuousProtonatableBead(id, "HCOOH", r, 1, spec, true);
            assignMomentum(bead, temperature);
        }
        
        std::clog << "Created " << cg->numberOfParticles() << " beads." << std::endl;
        std::clog << "Created " << cg->numberOfParticleGroups() << " waters groups" << std::endl;
        std::clog << "Created " << cg->numberOfFreeParticles() << " HCCOH beads." << std::endl;
        return cg;
    }
}