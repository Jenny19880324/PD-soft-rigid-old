// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SHAPE_MATRIX_H
#define IGL_SHAPE_MATRIX_H

#include "igl_inline.h"
#include "RBCEnergyType.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
 
namespace igl
{
  // Construct the shape matrix for a given rbc energy
  // Inputs:
  //   V  #V by Vdim list of initial domain positions
  //   F  #F by 4 list of triangle indices into V
  //   energy  RBCEnergyType enum value defining which energy is being used.
  //     See RBCEnergyType.h for valid options and explanations.
  // Outputs:
  //   SM dim*#F by dim sparse matrix, each 3x3 submatrix is the restpose shape matrix
  IGL_INLINE void shape_matrix(
    const Eigen::MatrixXd & V, 
    const Eigen::MatrixXi & F,
    const RBCEnergyType energy,
    Eigen::MatrixXd & SM);
}

#ifndef IGL_STATIC_LIBRARY
#include "shape_matrix.cpp"
#endif
#endif
