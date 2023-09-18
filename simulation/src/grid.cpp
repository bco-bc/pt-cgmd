/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 1 September 2022, 13:27h
 */

#include "simploce/simulation/grid.hpp"
#include "simploce/simulation/cell.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/util/logger.hpp"
#include <memory>
#include <map>
#include <cassert>
#include <utility>
#include <vector>

namespace simploce {
    namespace grid {

        /**
         * Creates cells.
         * @param spacing Requested cell spacing (edge length of each cell).
         * @param box Simulation box.
         * @return Cells, spacing in each direction (spacingX, spacingY, spacingZ), and number of cells in each direction (nx, ny, nz).
         */
        std::tuple<std::map<std::string, cell_ptr_t>,
                   real_t, real_t, real_t, int, int, int>
        setup_(const length_t& spacing, const box_ptr_t& box) {
            util::Logger logger("simploce::grid::setup_()");
            logger.trace("Entering.");

            logger.debug(std::to_string(spacing()) + ": Requested grid spacing (edge length of each cell).");

            // Box dimensions.
            auto lengthX = box->lengthX();
            auto lengthY = box->lengthY();
            auto lengthZ = box->lengthZ();

            // Number of cubes in each direction.
            auto nx = int(lengthX / spacing());
            auto ny = int(lengthY / spacing());
            auto nz = int(lengthZ / spacing());

            // Adjusted spacing in each direction.
            auto spacingX = lengthX / real_t(nx);
            auto spacingY = lengthY / real_t(ny);
            auto spacingZ = lengthZ / real_t(nz);

            // Create cells.
            std::map<std::string, cell_ptr_t> cells{};
            for (auto i = 0; i != nx; ++i) {
                for (auto j = 0; j != ny; ++j) {
                    for (auto k = 0; k != nz; ++k) {
                        Cell::location_t location{i,j,k};
                        auto x = (i + 0.5) * spacingX;
                        auto y = (j + 0.5) * spacingY;
                        auto z = (k + 0.5) * spacingZ;
                        position_t r{x, y, z};
                        auto cell = std::make_shared<Cell>(location, r);
                        auto key = cell->stringLocation();
                        auto pair = std::make_pair(key, cell);
                        cells.insert(pair);
                    }
                }
            }
            logger.debug(std::to_string(nx) + ": Number of cells in x-direction.");
            logger.debug(std::to_string(spacingX) + ": Cell spacing (edge length) in x-direction.");
            logger.debug(std::to_string(ny) + ": Number of cells in y-direction.");
            logger.debug(std::to_string(spacingY) + ": Cell spacing (edge length) in y-direction.");
            logger.debug(std::to_string(nz) + ": Number of cells in z-direction.");
            logger.debug(std::to_string(spacingZ) + ": Cell spacing (edge length) in z-direction.");
            logger.info(std::to_string(cells.size()) + ": Number of cells.");

            logger.trace("Leaving.");
            return std::make_tuple(cells,
                                   spacingX, spacingY, spacingZ,
                                   nx, ny, nz);
        }

        /**
         * Specifies neighboring cells based on a cutoff distance for cell pairs.
         * @param cells Cells.
         * @param cutoff Cutoff distance.
         * @param bc Boundary condition.
         * @return Neighbors.
         */
        std::map<std::string, std::vector<cell_ptr_t>>
        setupNeighbors_(const std::vector<cell_ptr_t>& cells,
                        const dist_t& cutoff,
                        const bc_ptr_t& bc) {
            util::Logger logger("simploce::pairlist::setupNeighbors_()");
            logger.trace("Entering.");

            logger.debug(std::to_string(cutoff()) + ": Cutoff distance for cell pairs.");
            logger.debug(std::to_string(cells.size()) + ": Number of cells.");

            // Neighboring cells for each cell in the grid.
            std::map<std::string, std::vector<cell_ptr_t>> neighbors;

            real_t rc2 = cutoff() * cutoff();
            logger.debug(std::to_string(rc2) + ": Square cutoff distance for cell pairs.");

            std::size_t ave{0};
            for (auto iter_i = cells.begin(); iter_i != cells.end(); ++iter_i) {
                const auto& central = *iter_i;
                auto ri = central->position();
                std::vector<cell_ptr_t> neighboring{};
                for (auto iter_j = cells.begin(); iter_j != cells.end(); ++iter_j) {
                    if (iter_j != iter_i) {
                        auto cell_j = *iter_j;
                        auto rj = cell_j->position();
                        auto rij = bc->apply(ri, rj);
                        auto rij2 = norm_square<real_t>(rij);
                        if (rij2 <= rc2) {
                            // std::clog << "Added." << std::endl;
                            neighboring.emplace_back(cell_j);
                        }
                    }
                }
                ave += neighboring.size();
                auto pair = std::make_pair(central->stringLocation(), neighboring);
                neighbors.insert(pair);
            }
            ave /= cells.size();
            logger.debug(std::to_string(ave) + ": Average number of neighboring cells.");

            logger.trace("Leaving.");
            return std::move(neighbors);
        }

        /**
         * Specifies neighboring cells for each cell in the grid. This assumes that the edge length of cells is
         * 0.5 * cutoff.
         * @param cells All cells.
         * @param nx Number of cells in x-direction.
         * @param ny Number of cells in y-direction.
         * @param nz Number of cells in z-direction.
         * @param limit Limit value for neighbor selection. The default value selects neighbors with cell index in
         * (i-limit, i+limit), where i refers to the current cell index in each direction.
         * @return Neighbors. Map with cell location (string) as key.
         */
        std::map<std::string, std::vector<cell_ptr_t>>
        setupNeighbors_(const std::map<std::string, cell_ptr_t>& cells,
                        std::size_t nx,
                        std::size_t ny,
                        std::size_t nz,
                        int limit = 2) {
            util::Logger logger("simploce::pairlist::setupNeighbors_()");
            logger.trace("Entering.");

            // Neighboring cells for each cell in the grid.
            std::map<std::string, std::vector<cell_ptr_t>> neighbors;

            std::size_t ave{0};
            for (const auto& c : cells) {
                std::vector<cell_ptr_t> neighboring{};
                auto central = c.second;
                auto location = central->location();
                auto i = int(std::get<0>(location));
                auto j = int(std::get<1>(location));
                auto k = int(std::get<2>(location));
                // Below employs simple periodic wrapping.
                for (auto ii = i - limit; ii <= i + limit; ++ii) {
                    auto iii = ii < 0 ? nx - 1 : ii;
                    iii = iii >= nx ? 0 : iii;
                    for (auto jj = j - limit; jj <= j + limit; ++jj) {
                        auto jjj = jj < 0 ? ny - 1 : jj;
                        jjj = jjj >= ny ? 0 : jjj;
                        for (auto kk = k - limit; kk <= k + limit; ++kk) {
                            auto kkk = kk < 0 ? nz - 1 : kk;
                            kkk = kkk >= nz ? 0 : kkk;
                            // Must exclude itself.
                            bool itself = (iii == i) && (jjj == j) && (kkk == k);
                            if (!itself) {
                                Cell::location_t loc{iii, jjj, kkk};
                                auto iter = cells.find(Cell::stringLocation(loc));
                                assert(iter != cells.end());
                                auto &neighbor = iter->second;
                                neighboring.emplace_back(neighbor);
                            }
                        }
                    }
                }
                ave += neighboring.size();
                auto pair = std::make_pair(central->stringLocation(), neighboring);
                neighbors.insert(pair);
            }
            ave /= cells.size();
            logger.info(std::to_string(ave) + ": Average number of neighboring cells.");

            logger.trace("Leaving.");
            return std::move(neighbors);
        }
    }

    Grid::Grid(box_ptr_t box, bc_ptr_t bc, dist_t cutoff) :
        box_{std::move(box)}, bc_{std::move(bc)}, cutoff_{cutoff} {

        auto spacing = 0.50 * cutoff_;
        auto setup = grid::setup_(spacing, box_);
        cells_ = std::move(std::get<0>(setup));
        spacingX_ = std::get<1>(setup);
        spacingY_ = std::get<2>(setup);
        spacingZ_ = std::get<3>(setup);
        spacing = spacingX_ > spacingY_ ? spacingX_ : spacingY_;
        spacing = spacing() > spacingZ_ ? spacing : spacingZ_;
        cutoff_ = 2.0 * spacing;
        nx_ = std::get<4>(setup);
        ny_ = std::get<5>(setup);
        nz_ = std::get<6>(setup);

        //auto cells = this->cells();
        //neighbors_ = grid::setupNeighbors_(cells, cutoffSR_, bc_);
        neighbors_ = grid::setupNeighbors_(cells_, nx_, ny_, nz_, 2);
    }

    void
    Grid::assignParticlesToCells(const simploce::bc_ptr_t &bc,
                                 const std::vector<p_ptr_t>& particles) {
        static util::Logger logger("simploce::Grid::assignParticlesToCells");
        logger.trace("Entering.");

        // Clear cells of particles.
        this->clear_();

        for (const auto& p : particles) {
            auto r = p->position();
            r = bc->placeInside(r);
            int i = std::floor(r[0] / spacingX_);
            int j = std::floor(r[1] / spacingY_);
            int k = std::floor(r[2] / spacingZ_);
            auto insideX = (i >= 0 && i < nx_);
            auto insideY = (j >= 0 && j < ny_);
            auto insideZ = (k >= 0 && k < nz_);
            //std::clog << p->index() << " " << i << " " << j << " " << k << " " << insideX << " " << insideY << " " << insideZ << std::endl;
            assert(insideX == true);
            assert(insideY == true);
            assert(insideZ == true);
            auto key = Cell::stringLocation(Cell::location_t{i,j,k});
            auto pair = cells_.find(key);
            auto cell = pair->second;
            cell->assign(p);
        }
        real_t ave{0};
        for (const auto& c : cells_) {
            auto& cell = c.second;
            ave += real_t(cell->particles().size());
        }
        ave /= real_t(cells_.size());
        logger.debug(std::to_string(particles.size()) +
                     ": Number of particles involved in non-bonded interactions.");
        logger.debug(std::to_string(ave) + ": Average number of particles in cells.");

        logger.trace("Leaving.");
    }

    const std::map<std::string, std::vector<cell_ptr_t>>&
    Grid::neighbors() const {
        return neighbors_;
    }

    cell_ptr_t
    Grid::findCell(const std::string& location) const {
        auto iter = cells_.find(location);
        assert(iter != cells_.end());
        return iter->second;
    }

    const std::vector<cell_ptr_t>&
    Grid::cells() const {
        static std::vector<cell_ptr_t> cells{};
        if (cells.empty()) {
            for (auto& iter: cells_) {
                auto cell = iter.second;
                cells.emplace_back(cell);
            }
        }
        return cells;
    }

    dist_t
    Grid::cutoff() const {
        return cutoff_;
    }

    void
    Grid::clear_() {
        for (auto& c : cells_) {
            auto& cell = c.second;
            cell->clear_();
        }
    }
}
