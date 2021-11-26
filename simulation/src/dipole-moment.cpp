/*
 * Author: André H. Juffer.
 * Created on 20/11/2021, 13:42.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/analysis/dipole-moment.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <string>

namespace simploce {

    DipoleMoment::DipoleMoment(stime_t dt,
                               real_t dm,
                               real_t maxStrengthGroup,
                               stime_t t0) :
            dt_{dt}, dm_{dm}, mmax_{maxStrengthGroup}, t0_{t0}, M_{}, hm_{} {
        auto nbins = std::size_t(maxStrengthGroup / dm);
        hm_.resize(nbins, 0);
    }

    void
    DipoleMoment::perform(const p_system_ptr_t& particleSystem) {
        static util::Logger logger("simploce::DipoleMoment::perform()");

        counter_ += 1;

        // Setup, if required.
        if (counter_ == 1) {
            particleSystem->doWithAllFreeGroups<void>([this](
                    const std::vector<p_ptr_t> &all,
                    const std::vector<p_ptr_t> &free,
                    const std::vector<pg_ptr_t> &groups) {
                auto maxStrength = this->mmax_;
                for (auto &g: groups) {
                    dipole_moment_t m{0.0, 0.0, 0.0};
                    for (auto &p: g->particles()) {
                        auto r = p->position();
                        auto q = p->charge();
                        for (std::size_t k = 0; k != 3; ++k) {
                            m[k] += q() * r[k];
                        }
                    }
                    auto strength = norm<real_t>(m);
                    if (strength > maxStrength) {
                        maxStrength = strength;
                    }
                }
                if (maxStrength > this->mmax_) {
                    this->mmax_ = 1.25 * maxStrength;
                    auto nbins = util::nint(real_t(this->mmax_) / this->dm_);
                    this->dm_ = mmax_ / real_t(nbins);
                    this->hm_.resize(nbins, 0);
                    std::string msg =
                            "Maximum value of -group- dipole moment strength m for "
                            "probability density function f(m) set to "
                            + util::toString(mmax_);
                    logger.debug(msg);
                }
                logger.debug("Bin size: " + util::toString(this->dm_));
            });
        }
        particleSystem->doWithAllFreeGroups<void>([this](
                const std::vector<p_ptr_t> &all,
                const std::vector<p_ptr_t> &free,
                const std::vector<pg_ptr_t> &groups) {
            dipole_moment_t M{0.0, 0.0, 0.0};
            for (const auto& p: all) {
                auto r = p->position();
                auto q = p->charge();
                for (std::size_t k = 0; k != 3; ++k) {
                    M[k] += q() * r[k];
                }
            }
            auto M2 = inner<real_t>(M, M);
            stime_t t = t0_ + counter_ * dt_;
            std::tuple<stime_t, dipole_moment_t, real_t> result = std::make_tuple(t, M, M2);
            this->M_.emplace_back(result);

            auto last = hm_.size() - 1;
            for (auto &g: groups) {
                dipole_moment_t m{0.0, 0.0, 0.0};
                for (auto &p: g->particles()) {
                    auto r = p->position();
                    auto q = p->charge();
                    for (std::size_t k = 0; k != 3; ++k) {
                        m[k] += q() * r[k];
                    }
                }
                auto strength = norm<real_t>(m);
                auto index = std::size_t (strength / dm_);
                if (index < (hm_.size() - 1)) {
                    hm_[index] += 1;
                } else {
                    hm_[last] += 1;
                }
            }
        });
    }

    DipoleMoment::result_t
    DipoleMoment::results() const
    {
        util::Logger logger("simploce::DipoleMoment::results");

        fm_result_t fm{};

        std::vector<real_t> hm(hm_.size(), 0.0);
        real_t tf = 0.0;
        for (std::size_t k = 0; k != hm_.size(); ++k) {
            hm[k] = real_t(hm_[k]) / real_t(counter_);  // Average over states.
            tf += real_t(hm_[k]);
        }
        std::size_t n = 0;
        for (auto v: hm) {
            real_t m = real_t(n) * dm_;
            real_t p = v / tf;          // Probability.
            real_t f = p / dm_;         // Probability density function.
            std::pair<real_t, real_t> pair = std::make_pair(m, f);
            fm.push_back(pair);
            n += 1;
        }

        logger.debug("Number of observations: " + util::toString(counter_));

        return std::make_pair(M_, fm);
    }

    dm_ptr_t
    DipoleMoment::create(const stime_t& dt,
                         real_t dm,
                         real_t maxGroupStrength,
                         const stime_t& t0) {
        return std::make_shared<DipoleMoment>(dt, dm, maxGroupStrength, t0);
    }

}