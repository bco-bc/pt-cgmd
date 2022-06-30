/*
 * File:   protonation-site-catalog.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 12, 2019, 5:35 PM
 */

#include "simploce/particle/protonation-site-catalog.hpp"
#include "simploce/particle/protonation-site.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/conf/p-conf.hpp"
#include <boost/algorithm/string.hpp>
#include <string>
#include <map>
#include <stdexcept>
#include <memory>

namespace simploce {
    
    // Atom iterator type.
    using atom_iter_t = std::vector<atom_ptr_t>::const_iterator;
    
    // Specifies a particle of a protonation site as listed in the catalog.
    struct catalog_particle_t {
        catalog_particle_t() : name{}, charge{0.0} {}
        
        std::string name;   // Particle name.
        charge_t charge;    // Charge value.
    };
    
    // Specifies a protonation site as listed in the catalog.
    struct catalog_protonation_site_t {
        catalog_protonation_site_t() : type{}, name{}, size{0}, deprotonated{}, protonated{} {}
        
        std::string type;                      // "type": Either "atom" or "bead"
        std::string name;                      // Particle name.
        std::size_t size;                      // Number of particles.
        std::vector<catalog_particle_t> deprotonated;  // Deprotonated state.
        std::vector<catalog_particle_t> protonated;    // protonated state.
    };
        
    static std::ostream& 
    operator << (std::ostream& stream, 
                 const catalog_protonation_site_t& site)
    {
        const int nameWidth = conf::NAME_WIDTH;
        
        stream << std::setw(nameWidth) << site.type << std::endl;
        stream << std::setw(nameWidth) << site.name << std::endl;
        stream << site.size << std::endl;
        for (std::size_t i = 0; i != site.deprotonated.size(); ++i) {
            stream << std::setw(nameWidth) << site.deprotonated[i].name;
            stream << std::setw(nameWidth) << site.deprotonated[i].charge;
            stream << std::setw(nameWidth) << site.protonated[i].charge;
            if ( i < (site.size - 1) ) {
                stream << std::endl;
            }
        }
        return stream;
    }
    
    /**
     * Selects atom for an actual protonation site.
     * @param first First identified atom. The requested atom should in the 
     * "neighbourhood" of the first atom.
     * @param window Search window around first atom.
     * @param name Name of atom sought.
     * @param atoms All atoms of physical system.
     * @return Iterator of selected atom in atoms, or atoms.end() if not found.
     */
    static atom_iter_t 
    selectAtom_(int first,
                int window,
                const std::string& name, 
                const std::vector<atom_ptr_t>& atoms)
    {        
        int natoms = int(atoms.size());
        
        // Define search window by assuming that all all atoms of protonation 
        // site are near the first identified atom (pointed to by <code>first</code>) 
        // of the site.
        int start = (first - window) >= 0 ? first - window : 0;       
        int end = (first + window) < natoms ? first + window : natoms;
        
        // Stop now?
        if ( start == end ) {
            return atoms.end();
        }
        
        // Find.
        int selected = start;
        bool found = false;
        int index = start;        
        while (!found && index <= end) {
            atom_ptr_t atom = atoms[index];
            if ( atom->name() == name ) {
                selected = index;
                found = true;
            }
            index += 1;
        }
        if ( found ) {           
           return atoms.begin() + selected;
        }
        if ( index == natoms ) {
            return atoms.end();  // Cannot go further.
        }
        
        // Try again with larger window.
        return selectAtom_(first, window + 5, name, atoms);
    }
    
    /**
     * Adds selected atom to protonation site, as currently specified by 
     * <code>particles</code>, <code>deprotonated</code>, and
     * <code>protonated</code>.
     * @param atom Previously selected atom.
     * @param index Index of particle in catalog protonation site.
     * @param site Catalog protonation site.
     * @param particles Protonation site particles.
     * @param deprotonated Deprotonated state protonation site.
     * @param protonated protonated state protonation site.
     */
    static void 
    addAtom_(const atom_ptr_t& atom,
             std::size_t index,
             const catalog_protonation_site_t& site,
             std::vector<atom_ptr_t>& particles,
             std::vector<spec_ptr_t>& deprotonated,
             std::vector<spec_ptr_t>& protonated)
    {
        // Add atom.
        particles.push_back(atom);
                        
        // Create new specification for deprotonated state.
        spec_ptr_t spec = atom->spec();
        bool free = spec->isFree();
        std::string dspecName = "D_" + spec->name();
        charge_t dcharge = site.deprotonated[index].charge;
        spec_ptr_t dspec = ParticleSpec::createFrom(spec, dspecName, dcharge, free, "Deprotonated.");
        deprotonated.push_back(dspec);  
        
        // Create new specification for protonated state.
        std::string pspecName = "P_" + spec->name();
        charge_t pcharge = site.protonated[index].charge;
        spec_ptr_t pspec = 
        ParticleSpec::createFrom(spec, pspecName, pcharge, free, "Protonated");
        protonated.push_back(pspec);                                
    }
    
    /*
     * Builds protonation sites according to catalog protonation site.
     * @param site Catalog protonation site.
     * @param atoms Atoms of physical system.
     * @return Protonation sites,
     */
    static std::vector<atom_prot_site_ptr_t> 
    buildAll_(const catalog_protonation_site_t& site, 
              const std::vector<atom_ptr_t>& atoms)
    {   
        // Initial search window.
        const int initial_window = 1;
        
        // All created protonation sites.
        std::vector<atom_prot_site_ptr_t> protonationSites{};
        
        // Name of first atom of given catalog protonation site.
        std::string nameFirstAtom = site.deprotonated[0].name;
        boost::trim(nameFirstAtom);

        std::size_t natoms = atoms.size();
        std::size_t counter = 0;
        bool found = false;
        
        // Protonation site.
        std::vector<atom_ptr_t> particles{};
        std::vector<spec_ptr_t> deprotonated{};
        std::vector<spec_ptr_t> protonated{};
                           
        // Build.
        while ( counter < natoms ) {
            
            // Identify next occurrence of first atom.
            do {
                atom_ptr_t atom = atoms[counter];
                found = (atom->name() == nameFirstAtom);
                counter += 1;
            } while ( !found && counter < (natoms - 1) );
                            
            if ( found ) {
                
                // Clear previous protonation site details.
                particles.clear();
                deprotonated.clear();
                protonated.clear();
                
                // Add first atom.
                auto first = counter - 1;
                atom_ptr_t atom = atoms[first];
                addAtom_(atom, 0, site, particles, deprotonated, protonated);
                                
                // Now find the rest of protonation site.
                for (std::size_t index = 1; index != site.size; ++index) {
                    std::string name = site.deprotonated[index].name;
                    boost::trim(name);
                    atom_iter_t selected = 
                            selectAtom_(first, initial_window, name, atoms);
                    if ( selected != atoms.end() ) {
                        atom = *selected;

                        // Add next atom.
                        addAtom_(atom, index, site, particles, deprotonated, protonated);                        
                    } else {
                        // This should never happen. If it happens, the provided
                        // atom collection of the physical system is faulty.
                        throw std::domain_error(
                            "Incomplete protonation site. "
                            "Cannot find the remaining constituting atoms."
                        );
                    }
                }

                std::vector<id_pair_t> bonds{};
                atom_prot_site_ptr_t protonationSite = 
                        std::make_shared<atom_prot_site_t>(site.name, 
                                                           particles,
                                                           bonds,
                                                           deprotonated,
                                                           protonated);
                protonationSites.push_back(protonationSite);
                
                // Continue.
                found = false;
            }                      
        }
        
        // Done.
        return protonationSites;
    }
    
    // Catalog
    using protonation_site_catalog_t = std::map<std::string, catalog_protonation_site_t>;
    
    static protonation_site_catalog_t catalog_{};
    
    
    std::vector<atom_prot_site_ptr_t> 
    ProtonationSiteCatalog::lookup(Atomistic& at) const
    {
        using prot_sites_cont_t = std::vector<atom_prot_site_ptr_t>;
        using iter_t = protonation_site_catalog_t::const_iterator;
        using pair_t = std::pair<std::string, catalog_protonation_site_t>;
        
        return at.doWithAll<prot_sites_cont_t>([] (const auto& atoms) {
            prot_sites_cont_t protonationSites{};
            
            for (iter_t iter = catalog_.begin(); iter != catalog_.end(); ++iter) {
                                
                // Get protonation site from catalog.
                pair_t pair = *iter;
                catalog_protonation_site_t site = pair.second;

                // Identify protonation sites.
                prot_sites_cont_t sites = buildAll_(site, atoms);
                for (auto ps : sites) {
                    protonationSites.push_back(ps);
                }
                                
            }
            return protonationSites;
        });
    }
    
    prot_site_catalog_ptr_t 
    ProtonationSiteCatalog::create(std::istream& stream)
    {
        using pair_t = std::pair<std::string, catalog_protonation_site_t>;
        
        const size_t bufferSize = conf::NAME_WIDTH;
        char charBuffer[bufferSize];
        std::string stringBuffer;
        
        // Skip first line.
        std::getline(stream, stringBuffer);  // Read EOL.
        
        std::string tname;
        std::getline(stream, tname);
        while ( stream.good() ) {
            catalog_protonation_site_t site;
            
            // Get the site type.
            boost::trim(tname);
            if ( tname != "atom" && tname != "bead") {
                throw std::domain_error(
                    tname + 
                    ": Illegal protonation site type. "
                    "Permissible values are 'atom' and 'bead'."
                );
            }
            site.type = tname;                        
            
            // Get the site name.            
            std::getline(stream, site.name);  // Read EOL.
            boost::trim(site.name);            
            
            // Get number of particles.
            std::size_t size;
            stream >> size;
            site.size = size;
            std::getline(stream, stringBuffer);  // Read EOL.
            
            real_t charge;
            catalog_particle_t protonated;
            catalog_particle_t deprotonated;
            for (std::size_t i = 0; i != size; ++i) {
                // Particle name
                stream.read(charBuffer, bufferSize);
                protonated.name = std::string(charBuffer, bufferSize);
                boost::trim(protonated.name);
                deprotonated.name = protonated.name;
                
                // Deprotonated charge
                stream >> charge;
                deprotonated.charge = charge;
                
                // Protonated charge.
                stream >> charge;
                protonated.charge = charge;
                
                site.deprotonated.push_back(deprotonated);
                site.protonated.push_back(protonated);
                
                std::getline(stream, stringBuffer);  // Read EOL.
            }
            pair_t pair = std::make_pair(site.name, site);
            catalog_.insert(pair);
            
            std::getline(stream, tname);  // Read EOL.
        }
        return prot_site_catalog_ptr_t(new ProtonationSiteCatalog);
    }
    
    ProtonationSiteCatalog::ProtonationSiteCatalog()
    {        
    }
    
    std::ostream& 
    operator << (std::ostream& stream, 
                 const ProtonationSiteCatalog& catalog)
    {
        using iter_t = protonation_site_catalog_t::const_iterator;
        using pair_t = std::pair<std::string, catalog_protonation_site_t>;
        std::size_t size = catalog_.size();
        std::size_t index = 0;
        for (iter_t iter = catalog_.begin(); iter != catalog_.end(); ++iter) {
            pair_t pair = *iter;
            auto site = pair.second;
            stream << site;
            if ( index < (size - 1) ) {
                stream << std::endl;
            }
            index += 1;
        }
        return stream;
    }
}

