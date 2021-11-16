/*
 * File:   grid-cells-test.cpp
 * Author: juffer
 *
 * Created on 7 November 2019, 21:46
 */

#include "simploce/simulation/cell.hpp"
#include "simploce/simulation/grid.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/util/util.hpp"
#include "simploce/simulation/s-properties.hpp"
#include <cstdlib>
#include <iostream>
#include <utility>
#include <cmath>

using namespace simploce;

using cell_t = Cell<Bead>;
using location_t = cell_t::location_t;
using grid_t = Grid<Bead>;


static void printLocation(const location_t& location)
{
    std::clog << "Cell location x: " << std::get<0>(location) << std::endl;
    std::clog << "Cell location y: " << std::get<1>(location) << std::endl;
    std::clog << "Cell location z: " << std::get<2>(location) << std::endl;    
}

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "grid-cells-test test 1" << std::endl;
    
    position_t r{0.5, 0.5, 0.5};
    location_t location{1,2,3};
    cell_t cell(r, location);

    printLocation(location);
    std::clog << "Cell position: " << cell.position() << std::endl;
    std::clog << "Number of free particles: " << cell.free().size() << std::endl;
    std::clog << "Number of particle groups: " << cell.groups().size() << std::endl;
    
    auto box = factory::cube(7.23);
    std::clog << "Box size: " << box->size() << std::endl;
    length_t sideLength = 0.5 * util::cutoffDistance(box);
    std::clog << "Requested side length: " << sideLength << std::endl;
    auto grid = grid_t::make(box, sideLength);
    std::clog << "Grid side length: " << grid->sideLength() << std::endl;
    auto ncells = grid->numberOfCells();
    std::clog << "Number of cells: " << ncells << std::endl;
    auto n1 = grid->numberOfCellPerDimension();
    std::clog << "Number of cell along one direction: " << n1 << std::endl;
    std::clog << "Cell positions and locations:" << std::endl;
    std::size_t counter = 0;
    for (auto& c : grid->cells()) {
        std::clog << counter;
        std::clog << " " << c.position();
        std::clog << " " << std::get<0>(c.location())
                  << " " << std::get<1>(c.location())
                  << " " << std::get<2>(c.location())
                  << std::endl;
        counter += 1;
    }
    
    location_t loc1{1, 1, 1};
    std::clog << "Neighbors of:" << std::endl;
    printLocation(loc1);
    auto neighbors = grid->neighbors(loc1);
    std::clog << "Number of neighbors: " << neighbors.size() << std::endl;
    for (auto& cell : neighbors) {
        std::clog << cell.position();
        std::clog << " " << std::get<0>(cell.location())
                  << " " << std::get<1>(cell.location())
                  << " " << std::get<2>(cell.location())
                  << std::endl;
    }
    
    location_t loc0{0, 0, 0};
    std::clog << "Neighbors of:" << std::endl;
    printLocation(loc0);
    neighbors = grid->neighbors(loc0);
    std::clog << "Number of neighbors: " << neighbors.size() << std::endl;
    for (auto& cell : neighbors) {
        std::clog << cell.position();
        std::clog << " " << std::get<0>(cell.location())
                  << " " << std::get<1>(cell.location())
                  << " " << std::get<2>(cell.location())
                  << std::endl;
    }
    
    std::clog << std::endl;
    auto& cell0 = grid->cell(loc0);
    auto r0 = cell0.position();
    std::clog << "Position cell0: " << r0 << std::endl;
    auto& cell1 = grid->cell(loc1);
    auto r1 = cell1.position();
    std::clog << "Position cell1: " << r1 << std::endl;

    auto bc = factory::pbc(box);
    auto rni = bc->apply(r0, r1);
    real_t RNI = norm<real_t>(rni);
    real_t R = norm<real_t>(r0 - r1);
    real_t rc = util::cutoffDistance(box)();
    std::clog << "Box size: " << box->size() << std::endl;
    std::clog << "Cutoff distance: " << rc << std::endl;
    std::clog << "Distance between (0,0,0) and (1,1,1): " << R << std::endl;
    std::clog << "Distance according to nearest image: " << RNI << std::endl;
    auto in = RNI < rc ? "YES" : "NO";
    std::clog << "Within cutoff distance?: " << in << std::endl;
}

void test2() {
    std::cout << "grid-cells-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (grid-cells-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% grid-cells-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (grid-cells-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (grid-cells-test)" << std::endl;

    //std::cout << "%TEST_STARTED% test2 (grid-cells-test)\n" << std::endl;
    //test2();
    //std::cout << "%TEST_FINISHED% time=0 test2 (grid-cells-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

