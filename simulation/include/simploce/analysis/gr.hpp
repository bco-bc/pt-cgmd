/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 22 September 2019, 15:22
 */

#ifndef GR_HPP
#define GR_HPP

#include "analyzer.hpp"
#include "a-types.hpp"
#include "simploce/types/s-types.hpp"
#include <string>
#include <memory>

namespace simploce {
    
    /**
     * Calculates the radial distribution function g(r) between particles of
     * given types.
     */
    class Gr : public analyzer {
    public:

        /**
         * Constructor
         * @param dr Spacing or bin size.
         * @param cutoff Cutoff distance for g(r) (maximal distance).
         * @param specName1 Particle specification name #1.
         * @param specName2 Particle specification name #2.
         * @param bc Boundary condition.
         */
        Gr(length_t dr,
           dist_t cutoff,
           std::string specName1,
           std::string specName2,
           bc_ptr_t bc);
        
        void 
        perform(const p_system_ptr_t& particleSystem) override;
    
        /**
         * Returns results.
         * @return Radial distribution function. First of each pair is r and the second
         * is g(r).
         */
        std::vector<std::pair<real_t, real_t>> 
        results() const;
        
        /**
         * Creates an analyzer.
         * @param dr Spacing or bin size.
         * @param specName1 Particle specification name #1.
         * @param specName2 Particle specification name #2.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return analyzer.
         */
        static gr_ptr_t
        create(const length_t& dr,
               const dist_t& cutoff,
               const std::string& specName1, 
               const std::string& specName2,
               const bc_ptr_t& bc);
    
    private:
        
        length_t dr_;
        dist_t cutoff_;
        std::string specName1_;
        std::string specName2_;
        bc_ptr_t bc_;
        
        // Helpers.
        void pp_(const std::vector<p_ptr_t>& particles);

        __attribute__((unused))
        void pg_(const std::vector<p_ptr_t>& particles,
                 const std::vector<pg_ptr_t>& groups);

        __attribute__((unused))
        void gg_(const std::vector<pg_ptr_t>& groups);
        
        std::size_t counter_{0};          // Counts number of states received.
        std::size_t nParticles1_{0};      // Number of particle of typeName #1.
        std::size_t nParticles2_{0};      // Number of particles of typeName #2.
        length_t rMax_{0.0};           // Upper limit g(r).
        volume_t volume_{0.0};         // Volume box.
        std::vector<std::size_t> hr_{};   // Counts number of particles in [r, r+dr].
    };
    

}

#endif /* GR_HPP */

