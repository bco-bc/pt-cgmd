/*
 * Author: André H. Juffer.
 * Created on 19/11/2021, 12:57.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/analysis/gr.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/math-constants.hpp"
#include <utility>
#include <cmath>

namespace simploce {

    Gr::Gr(length_t dr,
           dist_t cutoff,
           std::string specName1,
           std::string specName2,
           bc_ptr_t bc) :
        dr_{dr}, cutoff_{cutoff}, specName1_{std::move(specName1)}, specName2_{std::move(specName2)},
        bc_{std::move(bc)} {
        if ( specName1_.empty() || specName2_.empty() ) {
            throw std::domain_error("g(r): Two particle specification names must be provided.");
        }
        if ( !bc_ ) {
            throw std::domain_error("g(r): Missing boundary conditions.");
        }
    }

    void
    Gr::perform(const p_system_ptr_t& particleSystem)
    {
        static util::Logger logger("simploce::Gr::perform()");

        counter_ += 1;

        if ( counter_ == 1) {
            auto box = particleSystem->box();
            rMax_ = cutoff_;
            auto nBins = std::size_t(rMax_() / dr_());
            volume_ = box->volume();
            hr_.resize(nBins, 0);

            // Number of particles of the first and second specification.
            particleSystem->doWithAll<void>([this] (const std::vector<p_ptr_t>& all) {
                for (const auto& particle: all) {
                    if (particle->spec()->name() == this->specName1_) {
                        this->nParticles1_ += 1;
                    }
                    if (particle->spec()->name() == this->specName2_) {
                        this->nParticles2_ += 1;
                    }
                }
                this->pp_(all);

                logger.debug("Number of particles of typeName #1 '" + specName1_ + "': " + std::to_string(nParticles1_));
                logger.debug("Number of particles of typeName #2 '" + specName2_ + "': " + std::to_string(nParticles2_));
                logger.debug("Upper limit for g(r): " + std::to_string(rMax_()));
                logger.debug("Bin size: " + std::to_string(dr_()));
                if (nParticles1_ == 0 || nParticles2_ == 0) {
                    util::logAndThrow(logger, "No such particle(s).");
                }
            });
        } else {
            particleSystem->doWithAll<void>([this] (const std::vector<p_ptr_t>& all) {
                this->pp_(all);
            });
        }

      // Done.
    }

    void
    Gr::pp_(const std::vector<p_ptr_t>& particles)
    {
        static const real_t rc2 = rMax_() * rMax_();

        for (auto iter_i = particles.begin(); iter_i != particles.end(); ++iter_i) {
            const auto& pi = *iter_i;
            if ( pi->spec()->name() == specName1_ ) {
                auto ri = pi->position();
                for (auto iter_j = particles.begin(); iter_j != particles.end(); ++iter_j) {
                    if ( iter_i != iter_j) {
                        const auto& pj = *iter_j;
                        if ( pj->spec()->name() == specName2_ ) {
                            auto rj = pj->position();
                            auto rij = bc_->apply(ri, rj);
                            auto Rij2 = norm_square<real_t>(rij);
                            if ( Rij2 < rc2) {
                                auto Rij = std::sqrt(Rij2);
#if _DEBUG
                                util::tooClose<P>(pi, pj, std::tuple<energy_t, force_t, length_t>{});
#endif
                                auto index = std::size_t(Rij / dr_());
                                if ( index < hr_.size() ) {
                                    hr_[index] += 1.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    __attribute__((unused)) void
    Gr::pg_(const std::vector<p_ptr_t>& particles,
            const std::vector<pg_ptr_t>& groups) {
        for (const auto& pi : particles) {
            auto ri = pi->position();
            if ( pi->spec()->name() == specName1_ ) {
                for (const auto& g : groups) {
                    auto rg = g->position();
                    dist_vect_t r = bc_->apply(ri, rg);
                    length_t R = norm<length_t>(r);
                    if ( R() <= rMax_() ) {
                        for (const auto& pj : g->particles()) {
                            if ( pj->spec()->name() == specName2_ ) {
                                auto rj = pj->position();
                                dist_vect_t rij = bc_->apply(ri, rj);
                                length_t Rij = norm<length_t>(rij);
                                auto index = std::size_t(Rij() / dr_());
                                if ( index < hr_.size() ) {
                                    hr_[index] += 1.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    __attribute__((unused)) void
    Gr::gg_(const std::vector<pg_ptr_t>& groups)
    {
        for (auto it_i = groups.begin(); it_i != groups.end(); ++it_i) {
            const auto& gi = *it_i;
            auto particles_i = gi->particles();
            auto rgi = gi->position();
            for (const auto& gj : groups) {
                auto rgj = gj->position();
                auto rbc = bc_->apply(rgi, rgj);
                auto Rbc = norm<length_t>(rbc);
                if ( Rbc <= rMax_() ) {
                    auto particles_j = gj->particles();
                    for (const auto& particle_i : particles_i) {
                        if (particle_i->spec()->name() == specName1_) {
                            auto ri = particle_i->position();
                            for (const auto& particle_j : particles_j) {
                                if ( particle_j->spec()->name() == specName2_ &&
                                     particle_i != particle_j) {
                                    auto rj = particle_j->position();
                                    auto rij = bc_->apply(ri, rj);
                                    auto Rij = norm<length_t>(rij);
                                    auto index = std::size_t(Rij() / dr_());
                                    if ( index < hr_.size() ) {
                                        hr_[index] += 1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // See Friedman. "A course in statistical mechanics", Prentice Hall, 1985, p.82, Eq (4.20).
    std::vector<std::pair<real_t, real_t>>
    Gr::results() const {
        util::Logger logger{"simploce::Gr::results()"};

        std::vector<std::pair<real_t, real_t>> gr{};

        real_t factor = 4.0 * math::constants<real_t>::PI / 3.0;

        // Number density particle of the second specification.
        real_t rho2 = real_t(nParticles2_) / volume_();

        for ( std::size_t i = 0; i != hr_.size(); ++i) {

            // Compute volume dV of current shell.
            real_t ri = real_t(i) * dr_();
            real_t rii = real_t(( i + 1 )) * dr_();
            real_t dV = factor * ( rii * rii * rii - ri * ri * ri );

            // Number of particles of the second specification without correlation (ideal gas)
            real_t n_2 = rho2 * dV;

            // Normalise, and also average over the number particles of the
            // first specification and the number of states (observations).
            real_t g = counter_ > 0 && nParticles1_ > 0 ?
                       real_t(hr_[i]) / (n_2 * real_t(nParticles1_ * counter_)) :
                       0.0;

            // Save.
            auto pair = std::make_pair(ri, g);
            gr.push_back(pair);
        }

        logger.debug("Number of observations: " + std::to_string(counter_));

        return gr;
    }

    gr_ptr_t
    Gr::create(const length_t& dr,
               const dist_t& cutoff,
               const std::string& specName1,
               const std::string& specName2,
               const bc_ptr_t& bc) {
        return std::make_shared<Gr>(dr, cutoff, specName1, specName2, bc);
    }
}

