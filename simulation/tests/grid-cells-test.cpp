/*
 * The MIT License
 *
 * Copyright 2019 juffer.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   grid-cells-test.cpp
 * Author: juffer
 *
 * Created on 7 November 2019, 21:46
 */

#include "simploce/simulation/cell.hpp"
#include "simploce/simulation/grid.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/util/util.hpp"
#include "simploce/simulation/sim-util.hpp"
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
    for (auto& cell : grid->cells()) {
        std::clog << counter;
        std::clog << " " << cell.position();
        std::clog << " " << std::get<0>(cell.location())
                  << " " << std::get<1>(cell.location())
                  << " " << std::get<2>(cell.location())
                  << std::endl;
        counter += 1;
    }
    
    location_t loc2{1, 1, 1};
    std::clog << "Neighbors of:" << std::endl;
    printLocation(loc2);
    auto neighbors = grid->neighbors(loc2);
    std::clog << "Number of neighbors: " << neighbors.size() << std::endl;
    for (auto& cell : neighbors) {
        std::clog << cell.position();
        std::clog << " " << std::get<0>(cell.location())
                  << " " << std::get<1>(cell.location())
                  << " " << std::get<2>(cell.location())
                  << std::endl;
    }
    
    location_t loc3{0, 0, 0};
    std::clog << "Neighbors of:" << std::endl;
    printLocation(loc3);
    neighbors = grid->neighbors(location_t{0, 0, 0});
    std::clog << "Number of neighbors: " << neighbors.size() << std::endl;
    for (auto& cell : neighbors) {
        std::clog << cell.position();
        std::clog << " " << std::get<0>(cell.location())
                  << " " << std::get<1>(cell.location())
                  << " " << std::get<2>(cell.location())
                  << std::endl;
    }
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

