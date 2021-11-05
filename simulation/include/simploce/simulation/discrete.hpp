//
// Created by ajuffer on 10/28/21.
//

#ifndef SIMULATION_DISCRETE_HPP
#define SIMULATION_DISCRETE_HPP

#include "simploce/particle/protonatable.hpp"

namespace simploce {

    /**
     * Represents a discretely changing protonation state. Its values is 0 or 1.
     */
    class Discrete: public protonatable {
    public:

        Discrete();

        void protonate() override;

        void deprotonate() override;

        bool isProtonated() const;

        charge_t charge() const override;

        mass_t mass() const override;

    private:

        std::size_t x_;
    };
}

#endif //SIMULATION_DISCRETE_HPP
