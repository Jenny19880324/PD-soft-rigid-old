// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_FIT_RIGID_MOTION_H
#define IGL_FIT_RIGID_MOTION_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  // FIT_RIGID_MOTION computes the best-tting rigid transformation that aligns
  // two sets of corresponding points.
  // // Implementation from Least-Squares Rigid Motion Using SVD
  // Olga Sorkine-Hornung and Michael Rabinovich
  // https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
  //
  // Inputs:
  //   p  #V by dim, source vertices
  //   q  #V by dim, target vertices
  //   w  #V by 1, weight of each vertex
  // Outputs:
  //   R  dim by dim, the rotation matrix
  //   t  dim by 1, the translation vector
  template <typename DerivedV>
  IGL_INLINE void fit_rigid_motion(
    const Eigen::PlainObjectBase<DerivedV> & p,
    const Eigen::PlainObjectBase<DerivedV> & q,
	const Eigen::VectorXd & w,
          Eigen::Matrix3d & R,
		  Eigen::RowVector3d & t);
}

#ifndef IGL_STATIC_LIBRARY
#  include "fit_rigid_motion.cpp"
#endif

#endif
