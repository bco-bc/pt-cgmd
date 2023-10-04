//
// Created by ajuffer on 9/26/23.
//

#include "simploce/simulation/bc-util.hpp"
#include "simploce/util/direction.hpp"
#include "simploce/simulation/pbc-1d-sr.hpp"
#include <vector>

using namespace simploce;

int main() {

   Direction direction = Direction::valueOf('z');
   auto nc = bc::normalComponents(direction);
   std::cout << "Normal component indices for direction " << direction.value() << ": ";
   for (int k = 0; k != nc.size(); ++k) {
       std::cout << " " << nc[k];
   }
   std::cout << std::endl;

   box_ptr_t box = std::make_shared<box_t>(2.00000000e+01, 1.66666667e+01, 2.00000000e+01);
   position_t r1{0.1, 16.0, 15.0};
   std::cout << r1 << ", r1 crossed? " << bc::crossed(r1, box, direction) << std::endl;
   position_t r2{-1.0, 0, 0};
   std::cout << r2 << ", r2 crossed? " << bc::crossed(r2, box, direction) << std::endl;

   PBC_1D_SR bc{box, direction};
   velocity_t v1(1.0, 1.0, 1.0);
   std::cout << "Velocity before BC with r2: " << v1 << std::endl;
   auto v2 = bc.apply(v1, r2);
   std::cout << "Velocity after BC with r2: " << v2 << std::endl;

   std::cout << "Velocity before BC with r1: " << v1 << std::endl;
   v2 = bc.apply(v1, r1);
   std::cout << "Velocity after BC with r1: " << v2 << std::endl;

   position_t r3{21.0, 17.0, 15.0};
   std::cout << r3 << ", r3 crossed? " << bc::crossed(r3, box, direction) << std::endl;
   std::cout << "Velocity before BC with r3: " << v1 << std::endl;
   v2 = bc.apply(v1, r3);
   std::cout << "Velocity after BC with r3: " << v2 << std::endl;
}
