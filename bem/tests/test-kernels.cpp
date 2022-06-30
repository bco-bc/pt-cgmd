/*
 * Author: André H. Juffer.
 * Created on 24/05/2022, 23:54.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/kernels.hpp"
#include "simploce/bem/utils.hpp"
#include "simploce/types/cvector_t.hpp"
#include "simploce/util/math-constants.hpp"
#include "simploce/util/file.hpp"
#include <cmath>
#include <string>
#include <iostream>

using namespace simploce;

int main() {
    real_t eps1 = 20.0;
    real_t eps2 = 78.2;
    molarity_t I = 0.04;
    std::clog << "Ionic strength: " << I << " mol/l" << std::endl;
    auto ka = util::inverseDebijeLength<real_t>(I, eps2);
    std::clog << "Inverse Debije length: " << ka << " nm^-1" << std::endl;
    std::clog << "Debije length: " << 1.0/ka << " nm" << std::endl;
    real_t epsRatio = eps2 / eps1;
    position_t r{0.0, 0.0, 1.0};
    normal_t normal{0.0, 0.0, 1.0};
    position_t r0{r}; // Collocation point.
    normal_t normal0{r0[0], r0[1], r0[2]};

    std::string fn = "/wrk3/tests/kernels.dat";
    std::ofstream ostream;

    util::open_output_file(ostream, fn);
    int n = 100;
    real_t theta{0.0};
    for (std::size_t i = 1; i <= n; ++i) {
        theta += math::constants<real_t>::PI / real_t(n);
        r[1] = std::sin(theta);
        normal[1] = r[1];
        r[2] = std::cos(theta);
        normal[2] = r[2];
        auto lij0 = kernels::Lij0(epsRatio, r, normal, r0);
        auto kernels = kernels::KLMNij(ka, epsRatio, r,normal, r0, normal0);
        auto kij = std::get<0>(kernels);
        auto lij = std::get<1>(kernels);
        auto nij = std::get<2>(kernels);
        auto mij = std::get<3>(kernels);
        auto dis = norm<real_t>(r - r0);
        std::clog << "dis = " << dis << ", theta = " << theta << ", r,n(r) = " << r << " " << normal;
        std::clog << ", lij0, kij, lij, nij, mij = " << lij0 << " " << kij << " " << lij << " " << mij;
        std::clog << " " << nij << std::endl;
        ostream << dis << " " << lij0 << " " << kij << " " << lij << " " << mij << " " << nij << std::endl;
    }

    ostream.close();
}

