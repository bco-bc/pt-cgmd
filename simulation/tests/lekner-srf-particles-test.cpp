//
// Created by ajuffer on 2/6/24.
//

#include "simploce/potentials/lekner.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/util/file.hpp"
#include <iostream>
#include <fstream>

using namespace simploce;

int main() {
    auto box = factory::box(length_t{6.0}, length_t{6.0}, length_t{20.0});
    charge_t q2 = -0.01;
    auto negative =
            ParticleSpec::create("P-", q2, 1.0, 1.0, true, "NEG");
}