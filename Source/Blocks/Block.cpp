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
//#include "WaveAccumulationBlock.hpp"
#include "Block.hpp"
#elif defined(WITH_SOLVER_RUSANOV)
#include "Rusanov/RusanovBlock.hpp"
#endif
#endif

#include <cassert>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <omp.h>
#include <type_traits>
#include "WaveAccumulationBlock.hpp"
#include "WavePropagationBlock.hpp"

static constexpr RealType GRAVITY = 9.81f;

Blocks::Block1D::Block1D(const Tools::Float1D<RealType>& h_, const Tools::Float1D<RealType>& hu_, const Tools::Float1D<RealType>& hv_):
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

Blocks::Block* Blocks::Block::getBlockInstance(int nx, int ny, RealType dx, RealType dy, Tools::Float2D<RealType>& h, Tools::Float2D<RealType>& hu, Tools::Float2D<RealType>& hv) {
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

Blocks::Block::Block(int nx, int ny, RealType dx, RealType dy, Tools::Float2D<RealType>& h, Tools::Float2D<RealType>& hu, Tools::Float2D<RealType>& hv):
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
 // TODO THIS FUNCTION NEEDS TO PARALLELIZED
 // Initialize water height and discharge
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

 setBoundaryType(BoundaryEdge::Left, scenario.getBoundaryType(BoundaryEdge::Left));


 setBoundaryType(BoundaryEdge::Right, scenario.getBoundaryType(BoundaryEdge::Right));


 setBoundaryType(BoundaryEdge::Bottom, scenario.getBoundaryType(BoundaryEdge::Bottom));


 setBoundaryType(BoundaryEdge::Top, scenario.getBoundaryType(BoundaryEdge::Top));


 // Perform update after external write to variables

 synchAfterWrite();
}


void Blocks::Block::setWaterHeight(RealType (*h)(RealType, RealType)) {
 for (int i = 1; i <= nx_; i++) {
   for (int j = 1; j <= ny_; j++) {
     // Calculate the height of the cell using the h function
     h_[i][j] = h(offsetX_ + (i - RealType(0.5)) * dx_, offsetY_ + (j - RealType(0.5)) * dy_);
   }
 }


 synchWaterHeightAfterWrite();
}


void Blocks::Block::setDischarge(RealType (*u)(RealType, RealType), RealType (*v)(RealType, RealType)) {
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
 for (int i = 0; i <= nx_ + 1; i++) {
   for (int j = 0; j <= ny_ + 1; j++) {
     b_[i][j] = b;
   }
 }

 synchBathymetryAfterWrite();
}

void Blocks::Block::setBathymetry(RealType (*b)(RealType, RealType)) {
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
   //setGhostLayer();

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
 if (boundary_[BoundaryEdge::Left] == BoundaryType::Outflow || boundary_[BoundaryEdge::Left] == BoundaryType::Wall) {
   std::memcpy(b_[0], b_[1], sizeof(RealType) * (ny_ + 2));
 }

 if (boundary_[BoundaryEdge::Right] == BoundaryType::Outflow || boundary_[BoundaryEdge::Right] == BoundaryType::Wall) {
   std::memcpy(b_[nx_ + 1], b_[nx_], sizeof(RealType) * (ny_ + 2));
 }

 for (int i = 0; i <= nx_ + 1; i++) {

   if (boundary_[BoundaryEdge::Bottom] == BoundaryType::Outflow || boundary_[BoundaryEdge::Bottom] == BoundaryType::Wall) {
     b_[i][0] = b_[i][1];
   }
   if (boundary_[BoundaryEdge::Top] == BoundaryType::Outflow || boundary_[BoundaryEdge::Top] == BoundaryType::Wall) {
     b_[i][ny_ + 1] = b_[i][ny_];
   }
 }

 // Set corner values
 b_[0][0]             = b_[1][1];
 b_[0][ny_ + 1]       = b_[1][ny_];
 b_[nx_ + 1][0]       = b_[nx_][1];
 b_[nx_ + 1][ny_ + 1] = b_[nx_][ny_];

 // Synchronize after an external update of the bathymetry
 synchBathymetryAfterWrite();
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




void Blocks::Block::computeMaxTimeStep(const RealType dryTol, const RealType cfl) {
 // Initialize the maximum wave speed
 RealType maximumWaveSpeed = RealType(0.0);
 // Compute the maximum wave speed within the grid
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



// this is now useless
void Blocks::Block::initializeCornerGhostCells() {
 // Set values in corner ghost cells
 // ...
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


// Assuming this code is inside a method of some class

void Blocks::Block::applyBoundaryCondition(BoundaryEdge edge, int i) {
 bool negate;

 switch (edge) {
 case BoundaryEdge::Left:
   negate = (boundary_[edge] == BoundaryType::Wall);
   h_[0][i]  = h_[1][i];
   hu_[0][i] = (negate) ? -hu_[1][i] : hu_[1][i];
   hv_[0][i] = hv_[1][i];
   break;

 case BoundaryEdge::Right:
   negate = (boundary_[edge] == BoundaryType::Wall);
   h_[nx_ + 1][i]  = h_[nx_][i];
   hu_[nx_ + 1][i] = (negate) ? -hu_[nx_][i] : hu_[nx_][i];
   hv_[nx_ + 1][i] = hv_[nx_][i];
   break;

 case BoundaryEdge::Bottom:
   negate = (boundary_[edge] == BoundaryType::Wall);
   h_[i][0]  = h_[i][1];
   hu_[i][0] = hu_[i][1];
   hv_[i][0] = (negate) ? -hv_[i][1] : hv_[i][1];
   break;

 case BoundaryEdge::Top:
   negate = (boundary_[edge] == BoundaryType::Wall);
   h_[i][ny_ + 1]  = h_[i][ny_];
   hu_[i][ny_ + 1] = hu_[i][ny_];
   hv_[i][ny_ + 1] = (negate) ? -hv_[i][ny_] : hv_[i][ny_];
   break;
 default:
   // Handle default case or raise an error
   break;
 }
}

// this is now useless too
void Blocks::Block::setBoundaryConditions() {

 int end = ny_;
 for (int i = 0; i <= end; i++) {
   //left
   applyBoundaryCondition(BoundaryEdge::Left, i);
   // right
   applyBoundaryCondition(BoundaryEdge::Right, i);
 }

 end = nx_;

 for (int i = 0; i <= end; i++) {
   //bottom
   applyBoundaryCondition(BoundaryEdge::Bottom, i);
   // top
   applyBoundaryCondition(BoundaryEdge::Top, i);

 }

 initializeCornerGhostCells();
}

void Blocks::Block::applyBoundary(int x, int y, int i, BoundaryEdge edge) {
 h_[x][y]  = neighbour_[edge]->h[i];
 hu_[x][y] = neighbour_[edge]->hu[i];
 hv_[x][y] = neighbour_[edge]->hv[i];
}

// this is the function that calls all the other functions
void Blocks::Block::setGhostLayer() {
 // std::cout << "Set simple boundary conditions " << std::endl << std::flush;
 // Call to virtual function to set ghost layer values

 // merged into this function
 // setBoundaryConditions();

 printf("Hello im running");
 BoundaryType leftBoundary = boundary_[BoundaryEdge::Left];
 BoundaryType rightBoundary = boundary_[BoundaryEdge::Right];
 BoundaryType topBoundary = boundary_[BoundaryEdge::Top];
 BoundaryType bottomBoundary = boundary_[BoundaryEdge::Bottom];

 bool left = leftBoundary == BoundaryType::Connect;
 bool right = rightBoundary == BoundaryType::Connect;
 bool top = topBoundary == BoundaryType::Connect;
 bool bottom = bottomBoundary == BoundaryType::Connect;

 int end = std::max(nx_, ny_) + 1;

 bool initCorners = true;

#pragma omp parallel
#pragma omp single
#pragma omp taskloop grainsize(10)
 for (int i = 0; i <= end; i++) {
   if (i < ny_ && i > 0) {
     // Left
     applyBoundaryCondition(BoundaryEdge::Left, i);
     // Right
     applyBoundaryCondition(BoundaryEdge::Right, i);
   }

   if (i < nx_ && i > 0) {
     // Bottom
     applyBoundaryCondition(BoundaryEdge::Bottom, i);
     // Top
     applyBoundaryCondition(BoundaryEdge::Top, i);
   }

   if(initCorners){
     initializeCornerGhostCells();
     initCorners = false;
   }

   // setting left ghost layer
   if (left && i <= ny_ + 1 ) {
     applyBoundary(0, i, i, BoundaryEdge::Left);
   }
   // setting right ghost layer
   if (right && i <= ny_+ 1 ) {
     applyBoundary(nx_ + 1, i, i, BoundaryEdge::Right);
   }

   // setting bottom ghost layer
   if (bottom && i <= nx_+1 ) {
     applyBoundary(i, 0, i, BoundaryEdge::Bottom);
   }
   //  setting top ghost layer
   if (top && i <= nx_+1) {
     applyBoundary(i, ny_ + 1, i, BoundaryEdge::Top);
   }
 }


}


int Blocks::Block::getNx() const { return nx_; }

int Blocks::Block::getNy() const { return ny_; }
