// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_READJOINT_H
#define IGL_READJOINT_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace igl
{
  // load mesh joints from a .joint file
  //
  template <typename Index, typename DerivedV>
  IGL_INLINE bool readJOINT(
    const std::string joint_file_name,
	std::vector<std::vector<Index>> & I,
    Eigen::PlainObjectBase<DerivedV> & V);
  // Inputs:
  //   joint_file  pointer to already opened .joint file 
  // Outputs:
  //   mesh_file  closed file
  template <typename Index, typename DerivedV>
  IGL_INLINE bool readJOINT(
    FILE * joint_file,
	std::vector<std::vector<Index>> & I,
    Eigen::PlainObjectBase<DerivedV> & V);
}

#ifndef IGL_STATIC_LIBRARY
#  include "readMESH.cpp"
#endif

#endif
