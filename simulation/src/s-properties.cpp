/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 1:38 PM
 */

#include "simploce/simulation/s-properties.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/types/cvector_t.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/conf/s-conf.hpp"
#include <cmath>

namespace simploce {
    namespace properties {

        real_t kappa(const std::vector<p_ptr_t>& particles) {
            return 0.0;
        }

        temperature_t kineticTemperature(const std::vector<p_ptr_t>& particles,
                                         const energy_t& eKin,
                                         bool mesoscopic) {

            util::Logger logger("simploce::properties::kineticTemperature(...)");
            static int counter = 0;

            counter += 1;

            // Determine number of degrees of freedom (dof).
            std::size_t nParticles = 0;
            for (const auto& p: particles) {
                // Exclude frozen particles.
                nParticles += (p->frozen() ? 0 : 1);
            }
            real_t nDof = 3.0 * real_t(nParticles) - 3.0;  // Assuming total linear momentum is constant.
            if ( nDof < 3) {
                // No point in calculating the temperature for a low number of degrees of freedom.
                return 0.0;
            }

            // Subtract kinetic energy of center of mass (cm) motion.
            auto cmLinearMomentum = linearMomentum(particles);
            energy_t cmEKin, kineticEnergy = eKin;
            if (norm<real_t>(cmLinearMomentum) > 0.0) {
                mass_t totalMass = 0.0;
                for (auto& p: particles)
                    if ( !p->frozen() )
                        totalMass += p->mass();
                velocity_t cmVelocity{};
                for (size_t k = 0; k != 3; ++k) {
                    cmVelocity[k] = cmLinearMomentum[k] / totalMass();
                }
                cmEKin = 0.5 * totalMass() * inner<real_t>(cmVelocity, cmVelocity);
                kineticEnergy -= cmEKin;
            }

            if (counter == 1) {
                std::string s = "eKin: " + std::to_string(eKin()) +
                                ", cmEKin: " + std::to_string(cmEKin()) +
                                ", kineticEnergy: " + std::to_string(kineticEnergy());
                logger.debug(s);
                logger.debug(std::to_string(int(nDof)) + ": Number of degrees of freedom.");
            }

            // Ensure correct units.
            auto KB = mesoscopic ? 1 : units::mu<real_t>::KB;
            return 2.0 * kineticEnergy() / (nDof * KB);
       }

        pressure_t pressure(const std::vector<p_ptr_t>& particles,
                            const temperature_t& temperature,
                            const box_ptr_t& box,
                            bool mesoscopic)
        {
            volume_t volume = box->volume();
            real_t virial1 = 0.0;
            pressure_t pressure{};

            auto KB = mesoscopic ? 1 : units::mu<real_t>::KB;
            std::size_t nParticles = 0;
            for (const auto& particle : particles) {
                if (!particle->frozen()) {
                    position_t r = particle->position();
                    force_t f = particle->force();
                    virial1 += inner<real_t>(f, r);
                    nParticles += 1;
                }
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
                        "WARNING: Rij < " + std::to_string(rMin()) +
                        ", Rij = " + std::to_string(Rij()) +
                        ", pi = " + pi->name() + ", index = " + std::to_string(pi->index()) +
                        ", id = " + util::to_string(pi->id()) +
                        ", pj = " + pj->name() + ", index = " + std::to_string(pj->index()) +
                        ", id = " + util::to_string(pj->id()) +
                        ", energy: " + util::to_string(std::get<0>(efl));
                logger.warn(message);
            }
        }

                
    }
}