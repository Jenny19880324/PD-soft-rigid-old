// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_BONE_CONSTRAINT_TYPE_H
#define IGL_BONE_CONSTRAINT_TYPE_H
namespace igl
{
  enum BoneConstraintType
  {
    AFFINE_BONE_CONSTRAINT = 0,
	RIGID_BONE_CONSTRAINT = 1,
    NUM_BONE_CONSTRAINT_TYPES = 3
  };
}
#endif
