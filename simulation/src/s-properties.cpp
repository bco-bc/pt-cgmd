/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 1:38 PM
 */

#include "simploce/simulation/s-properties.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/types/cvector_t.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <cmath>
#include <limits>

namespace simploce {
    namespace properties {

        real_t kappa(const std::vector<p_ptr_t>& particles) {
            return 0.0;
        }

        temperature_t kineticTemperature(const std::vector<p_ptr_t>& particles,
                                         const energy_t& eKin,
                                         bool isMesoscale) {
            auto nParticles = particles.size();
            auto KB = isMesoscale ? 1 : units::mu<real_t>::KB;
            real_t nDof = 3.0 * real_t(nParticles) - 3.0;  // Assuming total linear momentum is constant.
            if (nDof > 3 ) {
                return 2.0 * eKin() / (nDof * KB);
            } else {
                // No point in calculating the temperature for a low number of degrees of freedom.
                return 0.0;
            }
        }

        pressure_t pressure(const std::vector<p_ptr_t>& particles,
                            const temperature_t& temperature,
                            const box_ptr_t& box,
                            bool isMesoscale)
        {
            volume_t volume = box->volume();
            real_t virial1 = 0.0;
            pressure_t pressure{};

            auto KB = isMesoscale ? 1 : units::mu<real_t>::KB;
            std::size_t nParticles = particles.size();
            for (const auto& particle : particles) {
                position_t r = particle->position();
                force_t f = particle->force();
                virial1 += inner<real_t>(f,r);
            }
            if ( volume() > 0.0 ) {
                virial1 /= ( 3.0 * volume() );
                real_t virial2 =
                        real_t(nParticles) * KB * temperature() / volume();
                pressure = virial2 - virial1; // In kJ/(mol nm^3)
            } else {
                pressure = 0.0;
            }
            return pressure;
        }

        real_t frohlich(real_t aveM2,
                        const temperature_t& temperature,
                        const box_ptr_t& box)
        {
            auto E0 = units::mu<real_t>::E0;
            auto kT = units::mu<real_t>::KB * temperature();
            auto volume = box->volume();
            real_t h = 1.0 / (E0 * volume) * aveM2 / (3.0 * kT);
            real_t a = 2.0;
            real_t b = -1.0 - 3.0 * h;
            real_t c = -1.0;
            real_t D = b * b - 4.0 * a * c;
            assert(D > 0);
            return (-b + std::sqrt(D)) / (2.0 * a);
        }

        static void
        tooClose(const p_ptr_t& pi,
                 const p_ptr_t& pj,
                 const std::tuple<energy_t, force_t, length_t>& efl)
        {
            static util::Logger logger("simploce::properties::tooClose");
            static length_t rMin = conf::SHORT;

            auto Rij = std::get<2>(efl);
            if (Rij() < rMin() ) {
                std::string message =
                        "WARNING: Rij < " + util::toString(rMin()) +
                        ", Rij = " + util::toString(Rij) +
                        ", pi = " + pi->name() + ", index = " + util::toString(pi->index()) +
                        ", id = " + util::toString(pi->id()) +
                        ", pj = " + pj->name() + ", index = " + util::toString(pj->index()) +
                        ", id = " + util::toString(pj->id()) +
                        ", energy: " + util::toString(std::get<0>(efl));
                logger.warn(message);
            }
        }

                
    }
}