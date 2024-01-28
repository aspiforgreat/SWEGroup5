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

// this is now useless
void Blocks::WavePropagationBlock::initializeCornerGhostCells() {
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

void Blocks::WavePropagationBlock::applyBoundaryCondition(BoundaryEdge edge, int i) {
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
void Blocks::WavePropagationBlock::setBoundaryConditions() {

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

void Blocks::WavePropagationBlock::applyBoundary(int x, int y, int i, BoundaryEdge edge) {
  h_[x][y]  = neighbour_[edge]->h[i];
  hu_[x][y] = neighbour_[edge]->hu[i];
  hv_[x][y] = neighbour_[edge]->hv[i];
}

// this is the function that calls all the other functions
void Blocks::WavePropagationBlock::setGhostLayer() {
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

void Blocks::WavePropagationBlock::computeVerticalEdgeUpdates(int i, int j, RealType& maxWaveSpeed) {
  if (j < ny_ + 1 && j > 0) {
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
}

void Blocks::WavePropagationBlock::computeHorizontalEdgeUpdates(int i, int j, RealType& maxWaveSpeed) {
  if (i < nx_ + 1 && i > 0) {
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

void Blocks::WavePropagationBlock::computeNumericalFluxes() {

  bool left = boundary_[BoundaryEdge::Left] == BoundaryType::Connect;
  bool right = boundary_[BoundaryEdge::Right] == BoundaryType::Connect;
  bool top = boundary_[BoundaryEdge::Top] == BoundaryType::Connect;
  bool bottom = boundary_[BoundaryEdge::Bottom] == BoundaryType::Connect;

  bool initCorners = true;
  bool settingBoundaries = true;

 // Maximum (linearized) wave speed within one iteration
 RealType maxWaveSpeed = RealType(0.0);
L
 //setGhostLayer();

 // Compute the net-updates for the vertical edges
 for (int i = 1; i < nx_ + 2; i++) {
    for (int j = 0; j < ny_ + 2; ++j) {

      if (settingBoundaries) {
        if (j < ny_ && j > 0) {
          // Left
          applyBoundaryCondition(BoundaryEdge::Left, j);
          // Right
          applyBoundaryCondition(BoundaryEdge::Right, j);
        }

        if (j < nx_ && j > 0) {
          // Bottom
          applyBoundaryCondition(BoundaryEdge::Bottom, j);
          // Top
          applyBoundaryCondition(BoundaryEdge::Top, j);
        }

        if (initCorners) {
          initializeCornerGhostCells();
          initCorners = false;
        }

        // setting left ghost layer
        if (left && j <= ny_ + 1) {
          applyBoundary(0, j, j, BoundaryEdge::Left);
        }
        // setting right ghost layer
        if (right && j <= ny_ + 1) {
          applyBoundary(nx_ + 1, j, j, BoundaryEdge::Right);
        }

        // setting bottom ghost layer
        if (bottom && j <= nx_ + 1) {
          applyBoundary(j, 0, j, BoundaryEdge::Bottom);
        }
        // setting top ghost layer
        if (top && j <= nx_ + 1) {
          applyBoundary(j, ny_ + 1, j, BoundaryEdge::Top);
        }
        settingBoundaries = false;
      }

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

    }// end of j loop

 } // end of i loop

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
