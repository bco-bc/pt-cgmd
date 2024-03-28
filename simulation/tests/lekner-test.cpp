/*
 * See also: A.H. Juffer, C.M. Shepherd, H.J. Vogel, "Protein–membrane electrostatic interactions:
 *           Application of the Lekner summation technique.", J. Chem. Phys., 114, 1892 - 1905,
 *           2001. http://dx.doi.org/10.1063/1.1334901
 *
 *  The test below reproduces Fig 2 of the reference above.
 *
 *  Created on 2 February 2024.
 */

#include "simploce/potentials/lekner.hpp"
#include "simploce/simulation/pbc-2d.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/util/file.hpp"
#include <iostream>
#include <fstream>

using namespace simploce;

int main() {

    real_t a = 0.2;
    auto box = factory::box(length_t{6.0}, length_t{6.0}, length_t{20.0});
    charge_t q1 = 0.4;
    charge_t q2 = -0.8;
    auto positive =
            ParticleSpec::create("P+", q1, 1.0, 1.0, true, "POS");
    auto negative =
            ParticleSpec::create("P-", q2, 1.0, 1.0, true, "NEG");
    auto test =
            ParticleSpec::create("T1", 1.0, 1.0, 1.0, true, "UT");

    Atomistic atomistic;

    // x- and y- coordinates of center of box parallel to XY-plane.
    auto cx = 0.5 * box->lengthX();
    auto cy = 0.5 * box->lengthY();

    // Three positive charges
    auto p1 = atomistic.addAtom("p1", positive);
    p1->position(position_t{cx-a, cy-a, 0.0});
    auto p2 = atomistic.addAtom("p2", positive);
    p2->position(position_t{cx+a, cy-a, 0.0});
    auto p3 = atomistic.addAtom("p3", positive);
    p3->position(position_t{cx-a, cy+a, 0.0});

    // One negative charge
    auto p4 = atomistic.addAtom("n1", negative);
    p4->position(position_t{cx+a, cy+a, 0.0});

    auto sigma = atomistic.charge() / (box->lengthX() * box->lengthY());
    std::cout << "Total reset: " << atomistic.charge() << std::endl;
    std::cout << "Surface reset density: " << sigma << " e/nm^2" << std::endl;

    // One test unit charge
    auto tp = Atom::create("a", 1, "TUC", test);
    tp->position(position_t{cx, cy, 0.0});

    std::string fileName = "/wrk3/tests/lekner.dat";
    std::ofstream ostream;
    util::open_output_file(ostream, fileName);

    bc_ptr_t bc = std::make_shared<PBC_2D>(box);
    Lekner lekner(box, bc);
    real_t dr = 0.01;
    int n = int(box->lengthZ() / (2.0 * dr));
    for (int i = 1; i < n; ++i) {
        real_t z = i * dr;
        tp->position(position_t{cx, cy, z});
        auto result1 = lekner.operator()(tp, p1);
        auto result2 = lekner.operator()(tp, p2);
        auto result3 = lekner.operator()(tp, p3);
        auto result4 = lekner.operator()(tp, p4);
        auto energy = result1.first + result2.first + result3.first + result4.first;
        std::cout << z << " " << energy << std::endl;
        ostream << z << " " << energy << std::endl;
    }
    ostream.close();
    std::cout << "Interaction energy written to \'" + fileName + "\'" << std::endl;
}