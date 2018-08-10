// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_RBCENERGYTYPE_H
#define IGL_RBCENERGYTYPE_H
namespace igl
{
  enum RBCEnergyType
  {
    RBC_ENERGY_TYPE_COROT = 0,
	RBC_ENERGY_TYPE_RBC = 1,
	RBC_ENERGY_TYPE_ACTIVE = 2, 
    RBC_ENERGY_TYPE_DEFAULT = 3,
    NUM_ARAP_ENERGY_TYPES = 4
  };
}
#endif
