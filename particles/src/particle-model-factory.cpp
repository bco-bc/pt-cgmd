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
#include <cmath>

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
        
    cg_ptr_t 
    ParticleModelFactory::harmonic(const length_t R0)
    {
        cg_ptr_t cg = std::make_shared<CoarseGrained>(); 
        auto spec = catalog_->lookup("AP");
        auto h = R0() / 2.0;
        position_t r1{0.0, 0.0, -h};
        auto bead_1 = cg->addBead(1, "AP1", r1, spec, false);
        position_t r2{0.0, 0.0, +h};
        auto bead_2 = cg->addBead(2, "AP2", r2, spec, false);
        std::vector<bead_ptr_t> beads{bead_1, bead_2};
        id_pair_t pair = std::make_pair(bead_1->id(), bead_2->id());
        std::vector<id_pair_t> pairs{pair};
        cg->addBeadGroup(beads, pairs);
        return cg;
    }
        
    cg_pol_water_ptr_t
    ParticleModelFactory::polarizableWater(const box_ptr_t& box,
                                           const density_t atDensitySI,
                                           const temperature_t temperature,
                                           bool protonatable,
                                           std::size_t nlimit) const
    {
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        std::clog << "Creating coarse grained particle model "
                     "for polarizable water." << std::endl;
        
        std::clog << "Temperature: " << temperature << std::endl;

        // Convert kg/m^3 to u/nm^3.
        density_t atDensity = atDensitySI / (SIUnits<real_t>::MU * 1.0e+27);
        std::clog << "Requested density (kg/m^3): " << atDensitySI << std::endl;
        std::clog << "Requested density (u/nm^3): " << atDensity << std::endl;
        
        // Spacing between DW water beads.
        const length_t spacing{0.53};  // Roughly the location of the first peak of g(r), in nm.
        std::clog << "Spacing between CW beads: " << spacing << " nm" << std::endl;
        
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
        std::size_t ncgWaters = real_t(natWaters) / 5.0;
        ncgWaters = ncgWaters > nlimit ? nlimit : ncgWaters;
        number_density_t cgNumberDensity = real_t(ncgWaters) / volume();
        std::size_t nOneDirection = std::pow(ncgWaters, 1.0/3.0);
        density_t cgDensity = cgNumberDensity * (cwSpec->mass() + dpSpec->mass());
        std::clog << "CG: Requested density (u/nm^3): " << cgDensity << std::endl;
        std::clog << "AT: Requested number density (1/nm^3): " << atNumberDensity << std::endl;
        std::clog << "AT: Requested number density (1/m^3): " << atNumberDensity*1.0e+27 << std::endl;
        std::clog << "CG: Requested number density (1/nm^3): " << cgNumberDensity << std::endl;
        std::clog << "CG: Requested number density (1/m^3): " << cgNumberDensity*1.0e+27 << std::endl;
        std::clog << "AT: Requested number of water molecules: " << natWaters << std::endl;
        std::clog << "CG: Requested number of water water groups: " << ncgWaters << std::endl;
        std::clog.flush();

        std::size_t nx = util::nint(Lx / spacing);
        std::size_t ny = util::nint(Ly / spacing);
        std::size_t nz = util::nint(Lz / spacing);
        std::clog << "Distance spacing between CW beads: " << spacing << std::endl;
        std::clog << "Number of coordinates in x-direction: " << nx << std::endl;
        std::clog << "Number of coordinates in y-direction: " << ny << std::endl;
        std::clog << "Number of coordinates in z-direction: " << nz << std::endl;
        std::clog << "Number of coordinates in any direction: " << nOneDirection << std::endl;
        std::clog.flush();
        
        // Start from an empty particle model.
        cg_pol_water_ptr_t cg = std::make_shared<PolarizableWater>();                
        
        // Add beads.
        real_t x0 = 0.0; //-0.5 * Lx();
        real_t y0 = 0.0; //-0.5 * Ly();
        real_t z0 = 0.0; //-0.5 * Lz();
        std::size_t counter = 0;
        std::size_t i = 0, j = 0, k = 0;         
        while ( i < nx && counter < ncgWaters ) {
            //real_t x = x0 + (i + 0.5) * spacing();
            real_t x = x0 + i * spacing();
            while ( j < ny && counter < ncgWaters ) {
                //real_t y = y0 + (j + 0.5) * spacing();
                real_t y = y0 + j * spacing();
                while ( k < nz && counter < ncgWaters ) {
                    //real_t z = z0 + (k + 0.5) * spacing();                    
                    real_t z = z0 + k * spacing();
                    std::size_t id = 2 * counter + 1;

                    // Single water group.
                    std::vector<bead_ptr_t> beads{};
                    std::vector<id_pair_t> bbonds{};
                            
                    
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
        if ( cg->numberOfParticleGroups() < ncgWaters ) {
            std::clog << "The number of created water groups is less then the number "
                         "of requested water groups." << std::endl;
        }
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
        std::size_t nHCOOH = util::nint(numberDensity() * volume());
        std::clog << "Number of HCOOH: " << nHCOOH << std::endl;
        
        std::clog << "Replacing " << nHCOOH << " water groups by HCOOH beads." << std::endl;
        spec_ptr_t spec = catalog_->lookup("HCOOH");
        auto size = cg->numberOfBeads();
        for (std::size_t counter = 0; counter != nHCOOH; ++counter) {
            std::size_t id = size + counter + 1;
            position_t r = cg->removeGroup_();
            auto bead = cg->addContinuousProtonatableBead(id, "HCOOH", r, 1, spec, true);
            if ( !protonatable ) {
                bead->deprotonate();
            }
            assignMomentum(bead, temperature);
        }
        
        std::clog << "Created " << cg->numberOfParticles() << " beads." << std::endl;
        std::clog << "Created " << cg->numberOfParticleGroups() << " waters groups" << std::endl;
        std::clog << "Created " << cg->numberOfFreeParticles() << " HCCOH beads." << std::endl;
        
        return cg;
    }
    
    cg_ptr_t 
    ParticleModelFactory::electrolyte(const box_ptr_t& box,
                                      molarity_t molarity,
                                      temperature_t temperature)
    {
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        std::clog << "Creating coarse grained particle model "
                     "for an electrolyte solution." << std::endl;
        
        std::clog << "Molarity (mol/l): " << molarity << std::endl;
        std::clog << "Temperature (K) : " << temperature << std::endl;


        spec_ptr_t NA = catalog_->lookup("Na+");
        spec_ptr_t CL = catalog_->lookup("Cl-");
        
        // Box details.
        length_t Lx = box->lengthX();
        length_t Ly = box->lengthY();
        length_t Lz = box->lengthZ();
        volume_t volume = box->volume();
        std::clog << "Unit Cell side length (nm): " << box->size() << std::endl;
        std::clog << "Unit Cell volume (nm^3): " << volume << std::endl;
        
        // Determine required number of ions.
        number_density_t numberDensity =
            molarity() * SIUnits<real_t>::NA / MUUnits<real_t>::l_to_nm3;
        std::size_t nions = util::nint(2.0 * numberDensity * volume);
        nions = (nions % 2 != 0 ? nions + 1 : nions);  // Neutrality condition.
        std::clog << "TOTAL number of ions requested: " << nions << std::endl;
    
        // Spacing between particles.
        const length_t spacing = 0.4; // The location of the first peak of g(r) for 
                                      // 0.1 M NaCl is around 0.3 nm.
        std::size_t nx = util::nint(Lx / spacing);
        std::size_t ny = util::nint(Ly / spacing);
        std::size_t nz = util::nint(Lz / spacing);
        std::clog << "Distance spacing between ion particles: " << spacing << std::endl;
        std::clog << "Number of coordinates in x-direction: " << nx << std::endl;
        std::clog << "Number of coordinates in y-direction: " << ny << std::endl;
        std::clog << "Number of coordinates in z-direction: " << nz << std::endl;
        
        // Start from an empty particle model.
        auto cg = std::make_shared<CoarseGrained>();
        
        // Create ions.
        std::size_t counter = 0;
        std::size_t i = 0, j = 0, k = 0;
        
        while ( i < nx && counter < nions ) {
            while ( j < ny && counter < nions ) {
                while ( k < nz && counter < nions ) {
                    // Id.
                    auto id = counter + 1;
                    
                    // Position.
                    real_t x = (i + 0.5) * spacing() + util::random<real_t>() * 0.1;
                    real_t y = (j + 0.5) * spacing() - util::random<real_t>() * 0.1;
                    real_t z = (k + 0.5) * spacing() + util::random<real_t>() * 0.1;
                    position_t r{x,y,z};
                    bead_ptr_t ion{};
                    if ( counter % 2 == 0 ) {
                        // Na+
                        ion = cg->addBead(id, NA->name(), r, NA, true);
                        assignMomentum<Bead>(ion, temperature);
                    } else {
                        // Na+
                        ion = cg->addBead(id, CL->name(), r, CL, true);
                        assignMomentum<Bead>(ion, temperature);
                    }                        
                    assignMomentum<Bead>(ion, temperature);
                    
                    counter += 1;
                    k += 1;
                }
                k = 0;
                j += 1;
            }
            j = 0;
            i += 1;
        }

        std::clog << "Number of ion particles generated: " << counter << std::endl;
        if ( counter < nions) {
            std::clog << "WARNING: Number of generated ions is smaller than requested.";
            std::clog << std::endl;
        }
        std::clog << "Ion number density (1/nm^3): "<< real_t(counter)/volume() << std::endl;
        
        return cg;        
    }
    
    cg_ptr_t 
    ParticleModelFactory::ljFluid(const box_ptr_t& box,
                                  const density_t atDensitySI,
                                  const temperature_t temperature)
    {
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        std::clog << "Creating coarse grained particle model "
                     "for LJ fluid." << std::endl;
        
        std::clog << "Temperature: " << temperature << std::endl;

        // Convert kg/m^3 to u/nm^3.
        density_t atDensity = atDensitySI / (SIUnits<real_t>::MU * 1.0e+27);
        std::clog << "Requested density (kg/m^3): " << atDensitySI << std::endl;
        std::clog << "Requested density (u/nm^3): " << atDensity << std::endl;
        
        // Spacing between LJ beads.
        const length_t spacing{0.5};
        std::clog << "Spacing between LJ beads: " << spacing << " nm" << std::endl;
        
        // Box details.
        length_t Lx = box->lengthX();
        length_t Ly = box->lengthY();
        length_t Lz = box->lengthZ();
        volume_t volume = box->volume();
        std::clog << "Box size (nm): " << box->size() << std::endl;
        std::clog << "Box volume (nm^3): " << volume << std::endl;
        
        std::size_t nbeads = atDensity() * volume();
        std::clog << "Requested number of LJ beads: " << nbeads << std::endl;
        
        std::size_t nx = util::nint(Lx / spacing);
        std::size_t ny = util::nint(Ly / spacing);
        std::size_t nz = util::nint(Lz / spacing);
        std::clog << "Number of coordinates in x-direction: " << nx << std::endl;
        std::clog << "Number of coordinates in y-direction: " << ny << std::endl;
        std::clog << "Number of coordinates in z-direction: " << nz << std::endl;
        
        // Start from an empty particle model.
        auto cg = std::make_shared<CoarseGrained>();
        
        // Apolar bead.
        auto spec = catalog_->lookup("AP");
        
        // Add beads.
        real_t x0 = 0.0;
        real_t y0 = 0.0;
        real_t z0 = 0.0;
        std::size_t counter = 0;
        std::size_t i = 0, j = 0, k = 0;         
        while ( i < nx && counter < nbeads ) {
            real_t x = x0 + (i + 0.5) * spacing();
            while ( j < ny && counter < nbeads ) {
                real_t y = y0 + (j + 0.5) * spacing();
                while ( k < nz && counter < nbeads ) {
                    real_t z = z0 + (k + 0.5) * spacing();                    
                    std::size_t id = counter + 1;
                    position_t r{x,y,z};
                    auto bead = cg->addBead(id, spec->name(), r, spec, true);
                    assignMomentum(bead, temperature);
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
        
        return cg;
    }
}