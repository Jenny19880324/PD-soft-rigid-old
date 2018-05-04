// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Alec Jacobson, Daniele Panozzo, Olga Diamanti 
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ROTATION_MATRIX_FROM_AXIS_AND_ANGLE_H
#define IGL_ROTATION_MATRIX_FROM_AXIS_AND_ANGLE_H
#include "igl_inline.h"

#include <Eigen/Core>

namespace igl 
{
  // Given rotation axis and angle calculate the rotation matrix from
  //
  // Inputs:
  //   u      3D column vector
  //   angle  Scalar
  // Output:
  //   3 by 3 rotation matrix that rotate angle around u axis
  //
  template <typename Scalar>
  IGL_INLINE Eigen::Matrix<Scalar, 3, 3> rotation_matrix_from_axis_and_angle(
    const Eigen::Matrix<Scalar, 3, 1> axis,
    const Scalar phi);
}

#ifndef IGL_STATIC_LIBRARY
#include "rotation_matrix_from_axis_and_angle.cpp"
#endif
#endif
