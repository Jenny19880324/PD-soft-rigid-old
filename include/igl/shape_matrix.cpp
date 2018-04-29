// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "shape_matrix.h"
#include "diag.h"
#include "sum.h"
#include "edges.h"
#include "verbose.h"
#include "cat.h"

IGL_INLINE void igl::shape_matrix(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const RBCEnergyType energy,
  Eigen::MatrixXd & SM)
{
  using namespace Eigen;
  // number of mesh vertices
  int n = V.rows();
  assert(n > F.maxCoeff());
  // dimension of mesh
  int dim = V.cols();
  // Number of mesh elements
  int m = F.rows();

  // number of rotations
  int nr;
  switch(energy)
  {
    case RBC_ENERGY_TYPE_COROT:
	case RBC_ENERGY_TYPE_RBC:
	case RBC_ENERGY_TYPE_DEFAULT:
      nr = m;
      break;
    default:
      fprintf(
        stderr,
        "shape_matrix.h: Error: Unsupported rbc energy %d\n",
        energy);
      return;
  }
  
  if(dim == 3)
  {
	  SM.resize(3 * m, dim);
	  for (int r = 0; r < m; r++)
	  {
		  Matrix3d shape_matrix;
		  for (int i = 0; i < 3; i++) {
			  for (int j = 0; j < 3; j++) {
				  shape_matrix(i, j) = V(F(r, i), j) - V(F(r, 3), j);
			  }
		  }
		  SM.block(r * 3, 0, 3, 3) = shape_matrix;
	  }
  }else
  {
    fprintf(
     stderr,
     "shape_matrix.h: Error: Unsupported dimension %d\n",
     dim);
    return;
  }

}
