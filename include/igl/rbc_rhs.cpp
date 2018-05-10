// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "rbc_rhs.h"
#include "verbose.h"
#include "repdiag.h"
#include "cat.h"
#include <iostream>

IGL_INLINE void igl::rbc_rhs(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const int dim,
  const igl::RBCEnergyType energy,
  Eigen::SparseMatrix<double>& J)
{
  using namespace std;
  using namespace Eigen;
  // Number of dimensions
  int Vdim = V.cols();

  switch(energy)
  {
    case RBC_ENERGY_TYPE_COROT:
	case RBC_ENERGY_TYPE_RBC:
      //nr = n;
      break;
    default:
      fprintf(
        stderr,
        "arap_rhs.h: Error: Unsupported rbc energy %d\n",
        energy);
      return;
  }

  J.resize(V.rows(), 3 * F.rows());
  J.setZero();
  std::vector<Eigen::Triplet<double>> J_triplets;
  J_triplets.reserve(12 * F.rows());
for(int i = 0; i < F.rows(); i++) {
	  	Matrix<double, 3, 3> Dm, Dm_inv_trans;
		Dm << V(F(i, 0), 0) - V(F(i, 3), 0), V(F(i, 1), 0) - V(F(i, 3), 0), V(F(i, 2), 0) - V(F(i, 3), 0),
			  V(F(i, 0), 1) - V(F(i, 3), 1), V(F(i, 1), 1) - V(F(i, 3), 1), V(F(i, 2), 1) - V(F(i, 3), 1),
			  V(F(i, 0), 2) - V(F(i, 3), 2), V(F(i, 1), 2) - V(F(i, 3), 2), V(F(i, 2), 2) - V(F(i, 3), 2);
		Dm_inv_trans = Dm.inverse().transpose();

		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				J_triplets.push_back(Eigen::Triplet<double>(F(i, row), 3 * i + col, Dm_inv_trans(col, row)));
			}
		}
		// row = 4
		for (int col = 0; col < 3; col++) {
			J_triplets.push_back(Eigen::Triplet<double>(F(i, 3), 3 * i + col,
				-Dm_inv_trans(col, 0) - Dm_inv_trans(col, 1) - Dm_inv_trans(col, 2)));
		}
  }

J.setFromTriplets(J_triplets.begin(), J_triplets.end());
J.makeCompressed();
}

