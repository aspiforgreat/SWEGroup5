/**
* @file
* This file is part of SWE.
*
* @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
* @author Sebastian Rettenberger (rettenbs AT in.tum.de,
* http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
* @author Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Michael_Bader)
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
* Implementation of Blocks::Block that uses solvers in the wave propagation formulation.
*/

#include "WavePropagationBlock.hpp"

#include <iostream>
#include <omp.h>

Blocks::WavePropagationBlock::WavePropagationBlock(int nx, int ny, RealType dx, RealType dy):
 Block(nx, ny, dx, dy),
 hNetUpdatesLeft_(nx + 1, ny),
 hNetUpdatesRight_(nx + 1, ny),
 huNetUpdatesLeft_(nx + 1, ny),
 huNetUpdatesRight_(nx + 1, ny),
 hNetUpdatesBelow_(nx, ny + 1),
 hNetUpdatesAbove_(nx, ny + 1),
 hvNetUpdatesBelow_(nx, ny + 1),
 hvNetUpdatesAbove_(nx, ny + 1) {}

Blocks::WavePropagationBlock::WavePropagationBlock(
 int nx, int ny, RealType dx, RealType dy,
 Tools::Float2D<RealType>& h,
 Tools::Float2D<RealType>& hu,
 Tools::Float2D<RealType>& hv
):
 Block(nx, ny, dx, dy, h, hu, hv),
 hNetUpdatesLeft_(nx + 1, ny),
 hNetUpdatesRight_(nx + 1, ny),
 huNetUpdatesLeft_(nx + 1, ny),
 huNetUpdatesRight_(nx + 1, ny),
 hNetUpdatesBelow_(nx, ny + 1),
 hNetUpdatesAbove_(nx, ny + 1),
 hvNetUpdatesBelow_(nx, ny + 1),
 hvNetUpdatesAbove_(nx, ny + 1) {}

void Blocks::WavePropagationBlock::setBoundary(const BoundaryEdge& edge, const std::function<void(bool, int, int, int, int)>& updateFunction) {
  int  start;
  int  end;
  bool leftRight = false;
  switch (edge) {
  case BoundaryEdge::Left:
  case BoundaryEdge::Right:
    leftRight = true;
    end       = ny_;
    break;
  case BoundaryEdge::Top:
  case BoundaryEdge::Bottom:
    // inccorect of course, just to understand
    leftRight = false;
    end       = nx_;
    break;
  }

  bool negate = false;
  switch (boundary_[edge]) {
  case BoundaryType::Wall:
    negate = true;
  case BoundaryType::Outflow:
    for (int j = 1; j <= end; j++) {

      if (leftRight) {
        updateFunction(negate, 0, j, nx_, ny_);
      } else {
        updateFunction(negate, j, 0, nx_, ny_);
      }
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

void Blocks::WavePropagationBlock::setLeftBoundary() {
  setBoundary(BoundaryEdge::Left, [&](bool negate, int i, int j, int nx_, int ny_) {
    h_[0][j]  = h_[1][j];
    hu_[0][j] = (negate) ? -hu_[1][j] : hu_[1][j];
    hv_[0][j] = hv_[1][j];
  });
}

void Blocks::WavePropagationBlock::setRightBoundary() {
  setBoundary(BoundaryEdge::Right, [&](bool negate, int i, int j, int nx_, int ny_) {
    h_[nx_ + 1][j]  = h_[nx_][j];
    hu_[nx_ + 1][j] = (negate) ? -hu_[nx_][j] : hu_[nx_][j];
    hv_[nx_ + 1][j] = hv_[nx_][j];
  });
}

void Blocks::WavePropagationBlock::setBottomBoundary() {
  setBoundary(BoundaryEdge::Bottom, [&](bool negate, int i, int j, int nx_, int ny_) {
    h_[i][0]  = h_[i][1];
    hu_[i][0] = hu_[i][1];
    hv_[i][0] = (negate) ? -hv_[i][1] : hv_[i][1];
  });
}

void Blocks::WavePropagationBlock::setTopBoundary() {
  setBoundary(BoundaryEdge::Top, [&](bool negate, int i, int j, int nx_, int ny_) {
    h_[i][ny_ + 1]  = h_[i][ny_];
    hu_[i][ny_ + 1] = hu_[i][ny_];
    hv_[i][ny_ + 1] = (negate) ? -hv_[i][ny_] : hv_[i][ny_];
  });
}

void Blocks::WavePropagationBlock::setBoundaryConditions() {
  // BoundaryType::Connect conditions are set in the calling function setGhostLayer
  // BoundaryType::Passive conditions need to be set by the component using Blocks::Block


  setLeftBoundary();


  setRightBoundary();


  setBottomBoundary();


  setTopBoundary();

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

// this is the function that calls all the other functions
void Blocks::WavePropagationBlock::setGhostLayer() {
  // std::cout << "Set simple boundary conditions " << std::endl << std::flush;
  // Call to virtual function to set ghost layer values

  setBoundaryConditions();

  // For a BoundaryType::Connect boundary, data will be copied from a neighbouring
  // Blocks::Block (via a Blocks::Block1D proxy object)
  // These copy operations cannot be executed in GPU/accelerator memory, e.g.,
  // setBoundaryConditions then has to take care that values are copied.

  // std::cout << "Set BoundaryType::Connect boundary conditions in main memory " << std::endl << std::flush;


  // Left boundary
  if (boundary_[BoundaryEdge::Left] == BoundaryType::Connect) {
    for (int j = 0; j <= ny_ + 1; j++) {
      h_[0][j]  = neighbour_[BoundaryEdge::Left]->h[j];
      hu_[0][j] = neighbour_[BoundaryEdge::Left]->hu[j];
      hv_[0][j] = neighbour_[BoundaryEdge::Left]->hv[j];
    }
  }

  // Right boundary
  if (boundary_[BoundaryEdge::Right] == BoundaryType::Connect) {

    for (int j = 0; j <= ny_ + 1; j++) {
      h_[nx_ + 1][j]  = neighbour_[BoundaryEdge::Right]->h[j];
      hu_[nx_ + 1][j] = neighbour_[BoundaryEdge::Right]->hu[j];
      hv_[nx_ + 1][j] = neighbour_[BoundaryEdge::Right]->hv[j];
    }
  }

  // Bottom boundary
  if (boundary_[BoundaryEdge::Bottom] == BoundaryType::Connect) {
    for (int i = 0; i <= nx_ + 1; i++) {
      h_[i][0]  = neighbour_[BoundaryEdge::Bottom]->h[i];
      hu_[i][0] = neighbour_[BoundaryEdge::Bottom]->hu[i];
      hv_[i][0] = neighbour_[BoundaryEdge::Bottom]->hv[i];
    }
  }

  // Top boundary
  if (boundary_[BoundaryEdge::Top] == BoundaryType::Connect) {

    for (int i = 0; i <= nx_ + 1; i++) {
      h_[i][ny_ + 1]  = neighbour_[BoundaryEdge::Top]->h[i];
      hu_[i][ny_ + 1] = neighbour_[BoundaryEdge::Top]->hu[i];
      hv_[i][ny_ + 1] = neighbour_[BoundaryEdge::Top]->hv[i];
    }
  }

}

void Blocks::WavePropagationBlock::computeNumericalFluxes() {
 // Maximum (linearized) wave speed within one iteration
 RealType maxWaveSpeed = RealType(0.0);

 // Compute the net-updates for the vertical edges
 for (int i = 1; i < nx_ + 2; i++) {
   for (int j = 1; j < ny_ + 2; ++j) {
     // printf("Hello from waveprog %d\n", omp_get_thread_num());
     if (j < ny_ + 1) {
       RealType maxEdgeSpeed = RealType(0.0);
       wavePropagationSolver_.computeNetUpdates(
         h_[i - 1][j],
         h_[i][j],
         hu_[i - 1][j],
         hu_[i][j],
         b_[i - 1][j],
         b_[i][j],
         hNetUpdatesLeft_[i - 1][j - 1],
         hNetUpdatesRight_[i - 1][j - 1],
         huNetUpdatesLeft_[i - 1][j - 1],
         huNetUpdatesRight_[i - 1][j - 1],
         maxEdgeSpeed
       );

       // Update the thread-local maximum wave speed
       maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);

     }

     // Compute the net-updates for the horizontal edges
     if (i < nx_ + 1) {
       RealType maxEdgeSpeed = RealType(0.0);
       wavePropagationSolver_.computeNetUpdates(
         h_[i][j - 1],
         h_[i][j],
         hv_[i][j - 1],
         hv_[i][j],
         b_[i][j - 1],
         b_[i][j],
         hNetUpdatesBelow_[i - 1][j - 1],
         hNetUpdatesAbove_[i - 1][j - 1],
         hvNetUpdatesBelow_[i - 1][j - 1],
         hvNetUpdatesAbove_[i - 1][j - 1],
         maxEdgeSpeed
       );

       // Update the thread-local maximum wave speed
       maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
     }
   }

 }

 if (maxWaveSpeed > 0.00001) {
   // Compute the time step width
   maxTimeStep_ = std::min(dx_ / maxWaveSpeed, dy_ / maxWaveSpeed);

   // Reduce maximum time step size by "safety factor"
   maxTimeStep_ *= RealType(0.4); // CFL-number = 0.5
 } else {
   // Might happen in dry cells
   maxTimeStep_ = std::numeric_limits<RealType>::max();
 }
}


void Blocks::WavePropagationBlock::updateUnknowns(RealType dt) {
 // Update cell averages with the net-updates
 for (int i = 1; i < nx_ + 1; i++) {
   for (int j = 1; j < ny_ + 1; j++) {
     h_[i][j] -= dt / dx_ * (hNetUpdatesRight_[i - 1][j - 1] + hNetUpdatesLeft_[i][j - 1])
                 + dt / dy_ * (hNetUpdatesAbove_[i - 1][j - 1] + hNetUpdatesBelow_[i - 1][j]);
     hu_[i][j] -= dt / dx_ * (huNetUpdatesRight_[i - 1][j - 1] + huNetUpdatesLeft_[i][j - 1]);
     hv_[i][j] -= dt / dy_ * (hvNetUpdatesAbove_[i - 1][j - 1] + hvNetUpdatesBelow_[i - 1][j]);

     if (h_[i][j] < 0) {
#ifndef NDEBUG
       // Only print this warning when debug is enabled
       // Otherwise we cannot vectorize this loop
       if (h_[i][j] < -0.1) {
         std::cerr << "Warning, negative height: (i,j)=(" << i << "," << j << ")=" << h_[i][j] << std::endl;
         std::cerr << "         b: " << b_[i][j] << std::endl;
       }
#endif

       // Zero (small) negative depths
       h_[i][j] = hu_[i][j] = hv_[i][j] = RealType(0.0);
     } else if (h_[i][j] < 0.1) {             // dryTol
       hu_[i][j] = hv_[i][j] = RealType(0.0); // No water, no speed!
     }
   }
 }
}
