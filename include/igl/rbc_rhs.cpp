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
for(int i = 0; i < F.rows(); i++) {
	  	Matrix<double, 3, 3> Dm, Dm_inv_trans;
		Dm << V(F(i, 0), 0) - V(F(i, 3), 0), V(F(i, 1), 0) - V(F(i, 3), 0), V(F(i, 2), 0) - V(F(i, 3), 0),
			  V(F(i, 0), 1) - V(F(i, 3), 1), V(F(i, 1), 1) - V(F(i, 3), 1), V(F(i, 2), 1) - V(F(i, 3), 1),
			  V(F(i, 0), 2) - V(F(i, 3), 2), V(F(i, 1), 2) - V(F(i, 3), 2), V(F(i, 2), 2) - V(F(i, 3), 2);
		Dm_inv_trans = Dm.inverse().transpose();
		
		SparseMatrix<double> Ai(3, V.rows());
		vector<Triplet<double>> Ai_IJV;
		Ai_IJV.push_back(Triplet<double>(0, F(i, 0), 1)); Ai_IJV.push_back(Triplet<double>(0, F(i, 3), -1));
		Ai_IJV.push_back(Triplet<double>(1, F(i, 1), 1)); Ai_IJV.push_back(Triplet<double>(1, F(i, 3), -1));
		Ai_IJV.push_back(Triplet<double>(2, F(i, 2), 1)); Ai_IJV.push_back(Triplet<double>(2, F(i, 3), -1));
		Ai.setFromTriplets(Ai_IJV.begin(), Ai_IJV.end());
		SparseMatrix<double> Gi(3, V.rows());
		Gi = (Dm_inv_trans * Ai).sparseView();
		
		SparseMatrix<double> Si(3, 3 * F.rows());
		vector<Triplet<double>> Si_IJV;
		Si_IJV.push_back(Triplet<double>(0, 3 * i + 0, 1));
		Si_IJV.push_back(Triplet<double>(1, 3 * i + 1, 1));
		Si_IJV.push_back(Triplet<double>(2, 3 * i + 2, 1));
		Si.setFromTriplets(Si_IJV.begin(), Si_IJV.end());
		J += Gi.transpose() * Si;
  }
}

