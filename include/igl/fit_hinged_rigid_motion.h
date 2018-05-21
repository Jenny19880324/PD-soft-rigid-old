// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_FIT_HINGED_RIGID_MOTION_H
#define IGL_FIT_HINGED_RIGID_MOTION_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
  // FIT_RIGID_MOTION computes the best-tting rigid transformation that aligns
  // two sets of corresponding points and satisfying a hinge constraint p.
  // // Implementation 
  // 1. drop the p2p constraint, find optimal rigid body motion using Procrustes
  // 2. linearize the constraint M(1)
  // 3. project M(1) using Procrustes
  // works two rigid bodies and one hinge now
  // Inputs:
  //   v1  #V1 of 1st rigid body by dim, restpose vertices
  //   d1  #V1 of 1st rigid body by dim, deformed vertices
  //   v2  #V2 of 2nd rigid body by dim, restpose vertices
  //   d2  #V2 of 2nd rigid body by dim, deformed vertices
  // Outputs:
  //   R1  dim by dim, the rotation matrix of the first rigid body
  //   t1  dim by 1, the translation vector of the first rigid body
  //   R2  dim by dim, the rotation matrix of the second rigid body
  //   t2  dim by 1, the translation vector of the second rigid body
  template <typename DerivedV>
  IGL_INLINE void fit_hinged_rigid_motion(
    const Eigen::PlainObjectBase<DerivedV> & v1,
	const Eigen::PlainObjectBase<DerivedV> & d1,
	const Eigen::PlainObjectBase<DerivedV> & cd1,
	const Eigen::PlainObjectBase<DerivedV> & v2,
	const Eigen::PlainObjectBase<DerivedV> & d2,
	const Eigen::PlainObjectBase<DerivedV> & cd2,
	const Eigen::RowVector3d &p,
    Eigen::Matrix3d &R1,
	Eigen::RowVector3d &t1,
	Eigen::Matrix3d &R2,
	Eigen::RowVector3d &t2);


  template <typename DerivedV>
  IGL_INLINE void fit_hinged_rigid_motion(
	  const Eigen::PlainObjectBase<DerivedV> & v1,
	  const Eigen::PlainObjectBase<DerivedV> & d1,
	  const Eigen::PlainObjectBase<DerivedV> & v2,
	  const Eigen::PlainObjectBase<DerivedV> & d2,
	  const Eigen::RowVector3d &p,
	  Eigen::Matrix3d &R1,
	  Eigen::RowVector3d &t1,
	  Eigen::Matrix3d &R2,
	  Eigen::RowVector3d &t2);


  template <typename DerivedV, typename DerivedN>
  IGL_INLINE void fit_hinged_rigid_motion(
	  const Eigen::PlainObjectBase<DerivedV> & V,
	  const Eigen::PlainObjectBase<DerivedN> & N,
	  const Eigen::PlainObjectBase<DerivedV> & P,
	  const std::vector<std::vector<int>> & I,
	  Eigen::PlainObjectBase<DerivedV> & U);
}

#ifndef IGL_STATIC_LIBRARY
#  include "fit_hinged_rigid_motion.cpp"
#endif

#endif
