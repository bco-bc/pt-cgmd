/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   model-factory.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 30, 2019, 12:51 PM
 */

#include "simploce/simulation/model-factory.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/cg-pol-water.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/util/mu-units.hpp"
#include "simploce/util/util.hpp"
#include <algorithm>
#include <random>
#include <cmath>
#include <memory>

namespace simploce {
    
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

    
    ModelFactory::ModelFactory(const spec_catalog_ptr_t& catalog) :
        catalog_{catalog}
    {        
    }
    
    
    cg_sim_model_ptr_t ModelFactory::createPolarizableWater(const box_ptr_t& box,
                                                            const density_t atDensitySI,
                                                            const temperature_t temperature)
    {
        std::clog.setf(std::ios_base::scientific, std::ios_base::floatfield);
        
        std::clog << "Creating Polarizable CGF Water:" << std::endl;

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
        
        // Periodic boundary conditions.
        bc_ptr_t bc = factory::pbc(box);
        std::clog << "Created periodic boundary conditions." << std::endl;
         
        // How many CG particles? A single polarizable CG water particle represents 
        // 5 water molecules. Each CG water particles consists of two connected 
        // CG particles (CW and DP).
        bead_spec_ptr_t mh2o = catalog_->molecularWater();
        bead_spec_ptr_t cw = catalog_->lookup("PCW");
        bead_spec_ptr_t dp = catalog_->lookup("DP");
        const length_t R_cw_dp = CoarseGrainedPolarizableWater::idealDistanceCWDP();
        std::clog << "\"Ideal\" distance between CW and DP: " << R_cw_dp 
                  << " nm" << std::endl;

        number_density_t atNumberDensity = atDensity / mh2o->mass();
        std::size_t natWaters = atNumberDensity * volume;
        std::size_t ncgWaters = util::nint(real_t(natWaters) / 5.0);
        number_density_t cgNumberDensity = real_t(ncgWaters) / volume();
        density_t cgDensity = cgNumberDensity * (cw->mass() + dp->mass());
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
    
        // Particle model.
        cg_ptr_t cg = std::make_shared<CoarseGrained>();
        
        // Add beads.
        std::size_t counter = 0;
        std::size_t i = 0, j = 0, k = 0;         
        while ( i < nx && counter < ncgWaters ) {
            while ( j < ny && counter < ncgWaters ) {
                while ( k < nz && counter < ncgWaters ) {

                    // Single water group.
                    std::vector<bead_ptr_t> beads{};
                    std::vector<id_pair_t> bbonds{};
                            
                    real_t x = (i + 0.5) * spacing();
                    real_t y = (j + 0.5) * spacing();
                    real_t z = (k + 0.5) * spacing();
                    position_t r1{x,y,z};
                    bead_ptr_t cwBead = cg->addBead("CW", r1, cw);
                    assignMomentum(cwBead, temperature);
                    beads.push_back(cwBead);
          
                    // Place the DP parallel to one of the 3 coordinate axes.
                    position_t r2{};
                    std::size_t l = util::random<real_t>() * 3.0;
                    switch (l) {
                        case 0: {
                            real_t coord = x + R_cw_dp();
                            r2 = position_t{coord, y, z};              
                            break;
                        }
                        case 1: {
                            real_t coord = y + R_cw_dp();
                            r2 = position_t{x,coord,z};
                            break;
                        }
                        default: {
                            real_t coord = z + R_cw_dp();
                            r2 = position_t{x,y,coord};
                            break;
                        }
                    } 
                    bead_ptr_t dpBead = cg->addBead("DP", r2, dp);
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
        std::clog << "Created " << cg->size() << " polarizable waters." << std::endl;

        // Interactor.
        cg_interactor_ptr_t interactor = 
                factory::interactorCoarseGrainedPolarizableWater(box, bc);
        
        // Displacer.
        std::shared_ptr<CoarseGrainedDisplacer> displacer = 
            std::make_shared<LangevinVelocityVerlet<CoarseGrained>>(interactor);
        
        // Done.
        return std::make_shared<cg_sim_model_t>(cg, displacer, interactor, box, bc);
    }    
    
    cg_sim_model_ptr_t ModelFactory::readCoarseGrainedFrom(std::istream& stream)
    {
        cg_sim_model_ptr_t cg{new SimulationModel<Bead>};
        cg->readFrom(stream, catalog_);
        return cg;
    }
  
}

