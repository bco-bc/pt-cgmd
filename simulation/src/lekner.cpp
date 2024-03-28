/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 2 February 2024. Original implementation copied from Reference 5 (listed in hardSphereLekner.hpp).
 */


#include "simploce/potentials/lekner.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/math-constants.hpp"
#include "simploce/util/bessel.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/conf/s-conf.hpp"
#include <cmath>

namespace simploce {
    namespace lekner {

        static real_t
        Bessel_K0(real_t x) {
            return math::Bessel_K0(x);
        }
    }

    Lekner::Lekner(simploce::box_ptr_t box,
                   bc_ptr_t bc,
                   real_t eps,
                   size_t n_max,
                   size_t k_max) :
       box_{std::move(box)}, bc_{std::move(bc)}, eps_{eps}, n_max_{n_max}, k_max_{k_max} {}

    std::pair<energy_t, force_t>
    Lekner::operator()(const simploce::p_ptr_t &pi, const simploce::p_ptr_t &pj) {
        static util::Logger logger("simploce::Lekner::operator()");
        logger.trace("Entering.");

        auto Rij = bc_->apply(pi->position(), pj->position());
        auto qi = pi->charge();
        auto qj = pj->charge();
        auto result = this->forceAndEnergy(Rij, qi, qj);

        logger.trace("Leaving.");
        return std::move(result);
    }

    std::pair<energy_t, force_t>
    Lekner::forceAndEnergy(const dist_vect_t& Rij,
                           const charge_t& qi,
                           const charge_t& qj) const {
        util::Logger logger("simploce::Lekner::forceAndEnergy()");static const real_t LARGE = conf::LARGE;
        logger.trace("Entering.");

        static auto two_pi = 2.0 * math::constants<real_t>::PI;
        static auto log_2 = std::log(2.0);
        static auto Lx = box_->lengthX();
        static auto Ly = box_->lengthY();
        static auto LyLyOverLxLx = Ly * Ly / (Lx * Lx);
        static auto LxLxOverLyLy = Lx * Lx / (Ly * Ly);

        auto R = norm<real_t>(Rij);
        auto dx = Rij[0];
        auto dy = Rij[1];
        auto dz = Rij[2];
        if (R < conf::SMALL) {
            logger.trace("Leaving.");
            return std::move(std::make_pair(energy_t{conf::LARGE}, force_t{}));
        }

        auto dxOverLx = dx / Lx;
        auto dyOverLy = dy / Ly;
        auto dzOverLy = dz / Ly;
        auto dzOverLx = dz / Lx;
        auto dzOverLx2 = dzOverLx * dzOverLx;
        auto dzOverLy2 = dzOverLy * dzOverLy;
        real_t energy;

        if (std::fabs(dy) > conf::SMALL) { // Whatever value of x and z.
            // Eq (7) of reference 5.
            // Sum over n.
            real_t d_sum_k;
            real_t d_sum_n;
            real_t sum_n = 0;
            int n = 1;
            do {
                auto sum_nP = sum_n;
                auto two_pi_n = two_pi * n;
                auto cos_term = std::cos(two_pi_n * dxOverLx);

                // Sum over k. Sum is initialized for k = 0.
                auto a2 = dyOverLy * dyOverLy;
                auto t1 = LyLyOverLxLx * a2 + dzOverLx2;
                auto t2 = two_pi_n * std::sqrt(t1);
                auto sum_k = lekner::Bessel_K0(t2);
                int k = 1;
                do {
                    auto sum_kP = sum_k;  // Current sum over k.
                    auto dk = real_t(k);

                    // Positive k.
                    auto a1 = dyOverLy + dk;
                    a2 = a1 * a1;
                    t1 = LyLyOverLxLx * a2 + dzOverLx2;
                    t2 = two_pi_n * std::sqrt(t1);
                    sum_k += lekner::Bessel_K0(t2);

                    // Negative k.
                    a1 = dyOverLy - dk;
                    a2 = a1 * a1;
                    t1 = LyLyOverLxLx * a2 + dzOverLx2;
                    t2 = two_pi_n * std::sqrt(t1);
                    sum_k += lekner::Bessel_K0(t2);

                    d_sum_k = sum_k - sum_kP;          // Difference, to test convergence.
                    k += 1;
                } while (std::fabs(d_sum_k) > this->eps_ && k < this->k_max_);
                sum_n += (cos_term * sum_k);
                d_sum_n = sum_n - sum_nP;             // Difference, to test convergence.
                n += 1;
            } while (std::fabs(d_sum_n) > this->eps_ && n < this->n_max_);

            // Convergence.
            auto qi_qj_over_Lx = qi * qj / Lx;
            energy = 4.0 * qi_qj_over_Lx * sum_n;
            auto a1 = std::cosh(two_pi * dzOverLy);
            auto a2 = std::cos(two_pi * dyOverLy);
            auto log_term = std::log(a1 - a2);
            energy -= (qi_qj_over_Lx * log_term);
            energy -= (qi_qj_over_Lx * log_2);    // Subtraction of just a constant.
            // Done.
        } else {                         // Whatever value of y and z.
            // Sum over k.
            // Eq (9) of reference 5.
            real_t sum_k = 0;
            real_t d_sum_n;
            real_t d_sum_k;
            int k = 1;
            do {
                auto sum_kP = sum_k;
                auto two_pi_k = two_pi * k;
                auto cos_term = std::cos(two_pi_k * dyOverLy);

                // Sum over n.  Sum is initialized for zero n.
                auto t1 = LxLxOverLyLy * dxOverLx * dxOverLx + dzOverLy2;
                auto t2 = two_pi_k*sqrt(t1);
                auto sum_n = lekner::Bessel_K0(t2);
                int n = 1;
                do {
                    auto sum_nP = sum_n;
                    auto dn = real_t(n);

                    // Positive n.
                    auto a1 = dxOverLx + dn;
                    auto a2 = a1 * a1;
                    t1 = LxLxOverLyLy * a2 + dzOverLy2;
                    t2 = two_pi_k * std::sqrt(t1);
                    sum_n += lekner::Bessel_K0(t2);

                    // Negative n.
                    a1 = dxOverLx - dn;
                    a2 = a1 * a1;
                    t1 = LxLxOverLyLy * a2 + dzOverLy2;
                    t2=two_pi_k*sqrt(t1);
                    sum_n += lekner::Bessel_K0(t2);

                    d_sum_n = sum_n - sum_nP;
                    n+=1;
                } while (std::fabs(d_sum_n) > this->eps_ && n < this->n_max_);

                sum_k += (cos_term * sum_n);
                d_sum_k = sum_k - sum_kP;
                k += 1;
            } while (std::fabs(d_sum_k) > this->eps_ && k < this->k_max_);

            // Convergence.
            auto qi_qj_over_Ly = qi * qj / Ly;
            energy = 4.0 * qi_qj_over_Ly * sum_k;
            auto a1 = std::cosh(two_pi * dzOverLx);
            auto a2 = std::cos(two_pi * dxOverLx);
            auto log_term = std::log(a1 - a2);
            energy -= (qi_qj_over_Ly * log_term);
            energy -= (qi_qj_over_Ly * log_2);     // Subtracting of just a constant.
        }
        energy /= units::mu<real_t>::FOUR_PI_E0;

        logger.trace("Leaving.");
        return std::move(std::make_pair(energy, force_t{0.0, 0.0, 0.0}));
    }

}