/*
 * File:   dipole-moment.hpp
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 24 October 2019, 21:21
 */

#ifndef DIPOLE_MOMENT_HPP
#define DIPOLE_MOMENT_HPP

#include "analyzer.hpp"
#include "a-types.hpp"
#include <utility>
#include <tuple>
#include <memory>

namespace simploce {
    
    /**
     * Computes the total dipole moment M and M*M at given times, as well as
     * the probability density function f(m) of the group dipole m.
     */
    class DipoleMoment : public Analyzer {
    public:

        /**
         * Holds result at a single time, time, M, and M*M.
         */
        using M_result_t = std::vector<std::tuple<stime_t, dipole_moment_t, real_t>>;
        
        /**
         * Holds the probability density f(m) of the strength or size of 
         * dipole moment m of particle groups. The pair is the pair m,f(m).
         */
        using fm_result_t = std::vector<std::pair<real_t, real_t>>;
        
        /**
         * Holds the result.
         */
        using result_t = std::pair<M_result_t, fm_result_t>;
        
        /**
         * Constructor
         * @param dt Time interval between successive states.
         * @param dm Bin size of probability density f(m).
         * @param maxStrengthGroup Maximum value of strength m (units are e nm) for f(m).
         * Default is 1.0 e nm.
         * @param t0 Start time. Default is 0.
         */
        DipoleMoment(stime_t dt,
                     real_t dm,
                     real_t maxStrengthGroup = 1.0,
                     stime_t t0 = 0);

        void perform(const p_system_ptr_t& particleSystem) override;
        
        /**
         * Returns results.
         * @return Results t (time), M and M*M, and m, f(m).
         */
        result_t 
        results() const;
        
        /**
         * Returns analyzer.
         * @param dt Time interval between successive states.
         * @param dm Bin size of probability density f(m).
         * @param maxGroupStrength Maximum value of strength m for f(m) (units are e nm). Default is
         * 10 e nm.
         * @param t0 Start time. Default is 0.
         * @return Analyzer.
         */
        static dm_ptr_t
        create(const stime_t& dt,
               real_t dm,
               real_t maxGroupStrength = 10,
               const stime_t& t0 = 0);
        
    private:
        
        stime_t dt_;    // Time step between states. Not the time step of the simulation.
        real_t dm_;     // Bin size
        real_t mmax_;   // Particle group dipole moment.
        stime_t t0_;    // Starting time.
        M_result_t M_;  // Total dipole moment.
        
        // Helpers.
        std::size_t counter_{0};
        std::vector<std::size_t> hm_;
    };
    
}

#endif /* DIPOLE_MOMENT_HPP */

