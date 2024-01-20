/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#include "Block.hpp"

#if !defined(ENABLE_CUDA)
#if defined(WITH_SOLVER_FWAVE) || defined(WITH_SOLVER_AUGRIE) || defined(WITH_SOLVER_HLLE)
#include "WaveAccumulationBlock.hpp"
#include "WavePropagationBlock.hpp"
#elif defined(WITH_SOLVER_RUSANOV)
#include "Rusanov/RusanovBlock.hpp"
#endif
#endif

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
#include <omp.h>
#include <functional>

static constexpr RealType GRAVITY = 9.81f;

Blocks::Block1D::Block1D(
  const Tools::Float1D<RealType>& h_, const Tools::Float1D<RealType>& hu_, const Tools::Float1D<RealType>& hv_
):
  h(h_),
  hu(hu_),
  hv(hv_){};

Blocks::Block1D::Block1D(RealType* h_, RealType* hu_, RealType* hv_, int size_, int stride_):
  h(h_, size_, stride_),
  hu(hu_, size_, stride_),
  hv(hv_, size_, stride_){};

#if defined(ENABLE_CUDA)
extern Block* getCUDABlockInstance(RealType, RealType, RealType, RealType);
#endif

Blocks::Block* Blocks::Block::getBlockInstance(int nx, int ny, RealType dx, RealType dy) {
  Block* block = nullptr;
#if !defined(ENABLE_CUDA)
#if defined(WITH_SOLVER_FWAVE) || defined(WITH_SOLVER_AUGRIE) || defined(WITH_SOLVER_HLLE)
  // block = new WaveAccumulationBlock(nx, ny, dx, dy);
  block = new WavePropagationBlock(nx, ny, dx, dy);
#elif defined(WITH_SOLVER_RUSANOV)
  block = new Rusanov::RusanovBlock(nx, ny, dx, dy);
#elif defined(WITH_SOLVER_AUGRIE_SIMD)
#error "Not implemented yet!"
#endif
#else
  block = getCUDABlockInstance(nx, ny, dx, dy);
#endif
  return block;
}

Blocks::Block* Blocks::Block::getBlockInstance(
  int nx, int ny, RealType dx, RealType dy,
  Tools::Float2D<RealType>& h,
  Tools::Float2D<RealType>& hu,
  Tools::Float2D<RealType>& hv) {
  Block* block = nullptr;
#if !defined(ENABLE_CUDA)
#if defined(WITH_SOLVER_FWAVE) || defined(WITH_SOLVER_AUGRIE) || defined(WITH_SOLVER_HLLE)
  // block = new WaveAccumulationBlock(nx, ny, dx, dy);
  block = new WavePropagationBlock(nx, ny, dx, dy, h, hu, hv);
#elif defined(WITH_SOLVER_RUSANOV)
  block = new Rusanov::RusanovBlock(nx, ny, dx, dy);
#elif defined(WITH_SOLVER_AUGRIE_SIMD)
#error "Not implemented yet!"
#endif
#else
  block = getCUDABlockInstance(nx, ny, dx, dy);
#endif
  return block;
}

Blocks::Block::Block(int nx, int ny, RealType dx, RealType dy):
  nx_(nx),
  ny_(ny),
  dx_(dx),
  dy_(dy),
  h_(nx + 2, ny + 2),
  hu_(nx + 2, ny + 2),
  hv_(nx + 2, ny + 2),
  b_(nx + 2, ny + 2),
  maxTimeStep_(0),
  offsetX_(0),
  offsetY_(0) {

  for (int i = 0; i < 4; i++) {
    boundary_[i]  = BoundaryType::Passive;
    neighbour_[i] = nullptr;
  }
}

Blocks::Block::Block(
  int nx, int ny, RealType dx, RealType dy,
  Tools::Float2D<RealType>& h,
  Tools::Float2D<RealType>& hu,
  Tools::Float2D<RealType>& hv):
  nx_(nx),
  ny_(ny),
  dx_(dx),
  dy_(dy),
  h_(h, true),
  hu_(hu, true),
  hv_(hv, true),
  b_(nx + 2, ny + 2),
  maxTimeStep_(0),
  offsetX_(0),
  offsetY_(0) {

  for (int i = 0; i < 4; i++) {
    boundary_[i]  = BoundaryType::Passive;
    neighbour_[i] = nullptr;
  }
}

void Blocks::Block::initialiseScenario(RealType offsetX, RealType offsetY, Scenarios::Scenario& scenario, const bool useMultipleBlocks) {
  offsetX_ = offsetX;
  offsetY_ = offsetY;

// Initialize water height and discharge
#pragma omp parallel
  {
    #pragma omp single
    {
      #pragma omp taskloop collapse(2)
      for (int i = 0; i <= nx_ + 1; i++) {
        for (int j = 0; j <= ny_ + 1; j++) {
          RealType x = offsetX + (i - RealType(0.5)) * dx_;
          RealType y = offsetY + (j - RealType(0.5)) * dy_;

          if (i >= 1 && i <= nx_ && j >= 1 && j <= ny_) {
            h_[i][j]  = scenario.getWaterHeight(x, y);
            hu_[i][j] = scenario.getVelocityU(x, y) * h_[i][j];
            hv_[i][j] = scenario.getVelocityV(x, y) * h_[i][j];
          }

          b_[i][j] = scenario.getBathymetry(x, y);
        }
      }

// Concurrently set boundary types
#pragma omp task
      setBoundaryType(BoundaryEdge::Left, scenario.getBoundaryType(BoundaryEdge::Left));

#pragma omp task
      setBoundaryType(BoundaryEdge::Right, scenario.getBoundaryType(BoundaryEdge::Right));

#pragma omp task
      setBoundaryType(BoundaryEdge::Bottom, scenario.getBoundaryType(BoundaryEdge::Bottom));

#pragma omp task
      setBoundaryType(BoundaryEdge::Top, scenario.getBoundaryType(BoundaryEdge::Top));
    }
  }

  // Perform update after external write to variables

  synchAfterWrite();
}


void Blocks::Block::setWaterHeight(RealType (*h)(RealType, RealType)) {
#pragma omp parallel
  {
#pragma omp single
    {
#pragma omp taskloop collapse(2)
      for (int i = 1; i <= nx_; i++) {
        for (int j = 1; j <= ny_; j++) {
          // Calculate the height of the cell using the h function
          h_[i][j] = h(offsetX_ + (i - RealType(0.5)) * dx_, offsetY_ + (j - RealType(0.5)) * dy_);
        }
      }
    }
#pragma omp taskwait
}


  synchWaterHeightAfterWrite();
}

// TODO find out if anyone uses this method
void Blocks::Block::setDischarge(RealType (*u)(RealType, RealType), RealType (*v)(RealType, RealType)) {
// TODO TASK
#pragma omp parallel for collapse(2) schedule(static)
  for (int i = 1; i <= nx_; i++) {
    for (int j = 1; j <= ny_; j++) {
      RealType x = offsetX_ + (i - RealType(0.5)) * dx_;
      RealType y = offsetY_ + (j - RealType(0.5)) * dy_;
      hu_[i][j]  = u(x, y) * h_[i][j];
      hv_[i][j]  = v(x, y) * h_[i][j];
    };
  }

  synchDischargeAfterWrite();
}

// TODO this method and the one below are not used anywhere, find out if they are needed
void Blocks::Block::setBathymetry(RealType b) {
  // TODO TASK
#pragma omp parallel for collapse(2) schedule(static)
  for (int i = 0; i <= nx_ + 1; i++) {
    for (int j = 0; j <= ny_ + 1; j++) {
      b_[i][j] = b;
    }
  }

  synchBathymetryAfterWrite();
}

void Blocks::Block::setBathymetry(RealType (*b)(RealType, RealType)) {
  // TODO TASK
#pragma omp parallel for collapse(2) schedule(static)
  for (int i = 0; i <= nx_ + 1; i++) {
    for (int j = 0; j <= ny_ + 1; j++) {
      b_[i][j] = b(offsetX_ + (i - RealType(0.5)) * dx_, offsetY_ + (j - RealType(0.5)) * dy_);
    }
  }

  synchBathymetryAfterWrite();
}

const Tools::Float2D<RealType>& Blocks::Block::getWaterHeight() {
  synchWaterHeightBeforeRead();
  return h_;
}

const Tools::Float2D<RealType>& Blocks::Block::getDischargeHu() {
  synchDischargeBeforeRead();
  return hu_;
}

const Tools::Float2D<RealType>& Blocks::Block::getDischargeHv() {
  synchDischargeBeforeRead();
  return hv_;
}

const Tools::Float2D<RealType>& Blocks::Block::getBathymetry() {
  synchBathymetryBeforeRead();
  return b_;
}

void Blocks::Block::simulateTimeStep(RealType dt) {
  computeNumericalFluxes();
  updateUnknowns(dt);
}

RealType Blocks::Block::simulate(RealType tStart, RealType tEnd) {
  RealType t = tStart;
  do {
    setGhostLayer();

    computeNumericalFluxes();
    updateUnknowns(maxTimeStep_);
    t += maxTimeStep_;

    std::cout << "Simulation at time " << t << std::endl << std::flush;
  } while (t < tEnd);

  return t;
}

void Blocks::Block::setBoundaryType(BoundaryEdge edge, BoundaryType boundaryType, const Block1D* inflow) {
  boundary_[edge]  = boundaryType;
  neighbour_[edge] = inflow;

  if (boundaryType == BoundaryType::Outflow || boundaryType == BoundaryType::Wall) {
    // One of the boundary was changed to BoundaryType::Outflow or BoundaryType::Wall
    // -> Update the bathymetry for this boundary
    setBoundaryBathymetry();
  }
}

void Blocks::Block::setBoundaryBathymetry() {
// Set bathymetry values in the ghost layer, if necessary
#pragma omp parallel
{
#pragma omp task
  {
    if (boundary_[BoundaryEdge::Left] == BoundaryType::Outflow || boundary_[BoundaryEdge::Left] == BoundaryType::Wall) {
      std::memcpy(b_[0], b_[1], sizeof(RealType) * (ny_ + 2));
    }
  }

#pragma omp task
  {
    if (boundary_[BoundaryEdge::Right] == BoundaryType::Outflow || boundary_[BoundaryEdge::Right] == BoundaryType::Wall) {
      std::memcpy(b_[nx_ + 1], b_[nx_], sizeof(RealType) * (ny_ + 2));
    }
  }
#pragma omp taskwait

#pragma omp for schedule(dynamic)
  for (int i = 0; i <= nx_ + 1; i++) {
    if (boundary_[BoundaryEdge::Bottom] == BoundaryType::Outflow || boundary_[BoundaryEdge::Bottom] == BoundaryType::Wall) {
      b_[i][0] = b_[i][1];
    }
    if (boundary_[BoundaryEdge::Top] == BoundaryType::Outflow || boundary_[BoundaryEdge::Top] == BoundaryType::Wall) {
      b_[i][ny_ + 1] = b_[i][ny_];
    }
  }

#pragma omp taskwait

  // Set corner values
  b_[0][0]             = b_[1][1];
  b_[0][ny_ + 1]       = b_[1][ny_];
  b_[nx_ + 1][0]       = b_[nx_][1];
  b_[nx_ + 1][ny_ + 1] = b_[nx_][ny_];

  // Synchronize after an external update of the bathymetry
  synchBathymetryAfterWrite();
}
}

  Blocks::Block1D* Blocks::Block::registerCopyLayer(BoundaryEdge edge) {
    switch (edge) {
    case BoundaryEdge::Left:
      return new Block1D(h_.getColProxy(1), hu_.getColProxy(1), hv_.getColProxy(1));
    case BoundaryEdge::Right:
      return new Block1D(h_.getColProxy(nx_), hu_.getColProxy(nx_), hv_.getColProxy(nx_));
    case BoundaryEdge::Bottom:
      return new Block1D(h_.getRowProxy(1), hu_.getRowProxy(1), hv_.getRowProxy(1));
    case BoundaryEdge::Top:
      return new Block1D(h_.getRowProxy(ny_), hu_.getRowProxy(ny_), hv_.getRowProxy(ny_));
    };
    return nullptr;
  }

Blocks::Block1D* Blocks::Block::grabGhostLayer(BoundaryEdge edge) {
  boundary_[edge] = BoundaryType::Passive;
  switch (edge) {
  case BoundaryEdge::Left:
    return new Block1D(h_.getColProxy(0), hu_.getColProxy(0), hv_.getColProxy(0));
  case BoundaryEdge::Right:
    return new Block1D(h_.getColProxy(nx_ + 1), hu_.getColProxy(nx_ + 1), hv_.getColProxy(nx_ + 1));
  case BoundaryEdge::Bottom:
    return new Block1D(h_.getRowProxy(0), hu_.getRowProxy(0), hv_.getRowProxy(0));
  case BoundaryEdge::Top:
    return new Block1D(h_.getRowProxy(ny_ + 1), hu_.getRowProxy(ny_ + 1), hv_.getRowProxy(ny_ + 1));
  };
  return nullptr;
}

void Blocks::Block::setGhostLayer() {
  // std::cout << "Set simple boundary conditions " << std::endl << std::flush;
  // Call to virtual function to set ghost layer values

  setBoundaryConditions();

  // For a BoundaryType::Connect boundary, data will be copied from a neighbouring
  // Blocks::Block (via a Blocks::Block1D proxy object)
  // These copy operations cannot be executed in GPU/accelerator memory, e.g.,
  // setBoundaryConditions then has to take care that values are copied.

  // std::cout << "Set BoundaryType::Connect boundary conditions in main memory " << std::endl << std::flush;


#pragma omp parallel
  {
#pragma omp single nowait
    {// Left boundary
     if (boundary_[BoundaryEdge::Left] == BoundaryType::Connect) {
#pragma omp task
       {
         for (int j = 0; j <= ny_ + 1; j++) {
           h_[0][j]  = neighbour_[BoundaryEdge::Left]->h[j];
          hu_[0][j] = neighbour_[BoundaryEdge::Left]->hu[j];
          hv_[0][j] = neighbour_[BoundaryEdge::Left]->hv[j];
        }
      }
    }

// Right boundary
if (boundary_[BoundaryEdge::Right] == BoundaryType::Connect) {
#pragma omp task
  {
    for (int j = 0; j <= ny_ + 1; j++) {
    h_[nx_ + 1][j]  = neighbour_[BoundaryEdge::Right]->h[j];
    hu_[nx_ + 1][j] = neighbour_[BoundaryEdge::Right]->hu[j];
    hv_[nx_ + 1][j] = neighbour_[BoundaryEdge::Right]->hv[j];
    }
  }
}

// Bottom boundary
if (boundary_[BoundaryEdge::Bottom] == BoundaryType::Connect) {
#pragma omp task
  {
    for (int i = 0; i <= nx_ + 1; i++) {
    h_[i][0]  = neighbour_[BoundaryEdge::Bottom]->h[i];
    hu_[i][0] = neighbour_[BoundaryEdge::Bottom]->hu[i];
    hv_[i][0] = neighbour_[BoundaryEdge::Bottom]->hv[i];
    }
  }
}

// Top boundary
if (boundary_[BoundaryEdge::Top] == BoundaryType::Connect) {
#pragma omp task
  {
    for (int i = 0; i <= nx_ + 1; i++) {
    h_[i][ny_ + 1]  = neighbour_[BoundaryEdge::Top]->h[i];
    hu_[i][ny_ + 1] = neighbour_[BoundaryEdge::Top]->hu[i];
    hv_[i][ny_ + 1] = neighbour_[BoundaryEdge::Top]->hv[i];
    }
  }
}
}
#pragma omp taskwait
}



            // std::cout << "Synchronize ghost layers (for heterogeneous memory) " << std::endl << std::flush;
  // Synchronize the ghost layers (for BoundaryType::Passive and BoundaryType::Connect conditions) with accelerator
  // memory
  synchGhostLayerAfterWrite();
}


void Blocks::Block::computeMaxTimeStep(const RealType dryTol, const RealType cfl) {
  // Initialize the maximum wave speed
  RealType maximumWaveSpeed = RealType(0.0);
  // Compute the maximum wave speed within the grid
#pragma omp parallel for collapse(2) reduction(max:maximumWaveSpeed) schedule(static)
  for (int i = 1; i <= nx_; i++) {
    for (int j = 1; j <= ny_; j++) {
      if (h_[i][j] > dryTol) {
        RealType momentum = std::max(std::abs(hu_[i][j]), std::abs(hv_[i][j]));

        RealType particleVelocity = momentum / h_[i][j];

        // Approximate the wave speed
        RealType waveSpeed = particleVelocity + std::sqrt(GRAVITY * h_[i][j]);

        maximumWaveSpeed = std::max(maximumWaveSpeed, waveSpeed);
      }
    }
  }

  RealType minimumCellLength = std::min(dx_, dy_);

  // Set the maximum time step variable
  maxTimeStep_ = minimumCellLength / maximumWaveSpeed;

  // Apply the CFL condition
  maxTimeStep_ *= cfl;
}

RealType Blocks::Block::getMaxTimeStep() const { return maxTimeStep_; }

void Blocks::Block::synchAfterWrite() {
  synchWaterHeightAfterWrite();
  synchDischargeAfterWrite();
  synchBathymetryAfterWrite();
}

void Blocks::Block::synchWaterHeightAfterWrite() {}

void Blocks::Block::synchDischargeAfterWrite() {}

void Blocks::Block::synchBathymetryAfterWrite() {}

void Blocks::Block::synchGhostLayerAfterWrite() {}

void Blocks::Block::synchBeforeRead() {
  synchWaterHeightBeforeRead();
  synchDischargeBeforeRead();
  synchBathymetryBeforeRead();
}

void Blocks::Block::synchWaterHeightBeforeRead() {}

void Blocks::Block::synchDischargeBeforeRead() {}

void Blocks::Block::synchBathymetryBeforeRead() {}

void Blocks::Block::synchCopyLayerBeforeRead() {}

void Blocks::Block::setBoundary(const BoundaryEdge& edge, const std::function<void(bool , int , int , int , int )>& updateFunction) {
  int start;
  int end;
  bool leftRight = false;
  switch (edge) {
    case BoundaryEdge::Left:
    case BoundaryEdge::Right:
        leftRight = true;
        end = ny_;
        break;
    case BoundaryEdge::Top:
    case BoundaryEdge::Bottom:
        //inccorect of course, just to understand
        leftRight = false;
        end = nx_;
        break;
  }

  bool negate = false;
  switch (boundary_[edge]) {
    case BoundaryType::Wall:
        negate = true;
    case BoundaryType::Outflow:
        #pragma omp parallel
          {
              #pragma omp single
              {
                    #pragma omp for schedule(dynamic)
                    for (int j = 1; j <= end; j++) {
                      #pragma omp task
                      {
                        if (leftRight) {
                          updateFunction(negate, 0, j, nx_, ny_);
                        } else {
                          updateFunction(negate, j, 0, nx_, ny_);
                        }
                      }
                    }
              }
              #pragma omp taskwait
          }
    break;


    case BoundaryType::Connect:
    case BoundaryType::Passive:
      break;
    default:
      assert(false);
      break;
  }

}

void Blocks::Block::setLeftBoundary() {
  setBoundary(BoundaryEdge::Left, [&](bool negate, int i, int j, int nx_, int ny_) {
    h_[0][j]  = h_[1][j];
    hu_[0][j] = (negate) ?  -hu_[1][j] : hu_[1][j] ;
    hv_[0][j] = hv_[1][j];
  });
}

void Blocks::Block::setRightBoundary() {
  setBoundary(BoundaryEdge::Right, [&](bool negate, int i, int j, int nx_, int ny_) {
    h_[nx_ + 1][j] = h_[nx_][j];
    hu_[nx_ + 1][j] = (negate) ? -hu_[nx_][j] : hu_[nx_][j];
    hv_[nx_ + 1][j] = hv_[nx_][j];
  });
}

void Blocks::Block::setBottomBoundary() {
  setBoundary(BoundaryEdge::Bottom, [&](bool negate, int i, int j, int nx_, int ny_) {
    h_[i][0]  = h_[i][1];
    hu_[i][0] = hu_[i][1];
    hv_[i][0] = (negate) ? -hv_[i][1] :  hv_[i][1] ;
  });
}

void Blocks::Block::setTopBoundary() {
  setBoundary(BoundaryEdge::Top, [&](bool negate, int i, int j, int nx_, int ny_ ) {
    h_[i][ny_ + 1]  = h_[i][ny_];
    hu_[i][ny_ + 1] = hu_[i][ny_];
    hv_[i][ny_ + 1] = (negate) ? -hv_[i][ny_] : hv_[i][ny_];
  });
}

void Blocks::Block::setBoundaryConditions() {
  // BoundaryType::Connect conditions are set in the calling function setGhostLayer
  // BoundaryType::Passive conditions need to be set by the component using Blocks::Block
#pragma omp parallel
  {
      #pragma omp single
      {
        #pragma omp task
                    setLeftBoundary();

          #pragma omp task
                    setRightBoundary();

          #pragma omp task
                    setBottomBoundary();

          #pragma omp task
                    setTopBoundary();

#pragma omp taskwait
      }
  }
  /*
   * Set values in corner ghost cells. Required for dimensional splitting and visualization.
   *   The quantities in the corner ghost cells are chosen to generate a zero Riemann solutions
   *   (steady state) with the neighboring cells. For the lower left corner (0,0) using
   *   the values of (1,1) generates a steady state (zero) Riemann problem for (0,0) - (0,1) and
   *   (0,0) - (1,0) for both outflow and reflecting boundary conditions.
   *
   *   Remark: Unsplit methods don't need corner values.
   *
   * Sketch (reflecting boundary conditions, lower left corner):
   * <pre>
   *                  **************************
   *                  *  _    _    *  _    _   *
   *  Ghost           * |  h   |   * |  h   |  *
   *  cell    ------> * | -hu  |   * |  hu  |  * <------ Cell (1,1) inside the domain
   *  (0,1)           * |_ hv _|   * |_ hv _|  *
   *                  *            *           *
   *                  **************************
   *                  *  _    _    *  _    _   *
   *   Corner Ghost   * |  h   |   * |  h   |  *
   *   cell   ------> * |  hu  |   * |  hu  |  * <----- Ghost cell (1,0)
   *   (0,0)          * |_ hv _|   * |_-hv _|  *
   *                  *            *           *
   *                  **************************
   * </pre>
   */
  h_[0][0]  = h_[1][1];
  hu_[0][0] = hu_[1][1];
  hv_[0][0] = hv_[1][1];

  h_[0][ny_ + 1]  = h_[1][ny_];
  hu_[0][ny_ + 1] = hu_[1][ny_];
  hv_[0][ny_ + 1] = hv_[1][ny_];

  h_[nx_ + 1][0]  = h_[nx_][1];
  hu_[nx_ + 1][0] = hu_[nx_][1];
  hv_[nx_ + 1][0] = hv_[nx_][1];

  h_[nx_ + 1][ny_ + 1]  = h_[nx_][ny_];
  hu_[nx_ + 1][ny_ + 1] = hu_[nx_][ny_];
  hv_[nx_ + 1][ny_ + 1] = hv_[nx_][ny_];
}

int Blocks::Block::getNx() const { return nx_; }

int Blocks::Block::getNy() const { return ny_; }