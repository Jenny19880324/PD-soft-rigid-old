// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "laplacian_matrix.h"
#include <vector>

// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>


template <typename DerivedV, typename DerivedF, typename Scalar>
IGL_INLINE void igl::laplacian_matrix(
	const Eigen::MatrixBase<DerivedV> & V,
	const Eigen::MatrixBase<DerivedF> & F,
	Eigen::SparseMatrix<Scalar>& L)
{
	using namespace Eigen;
	using namespace std;

	L.resize(V.rows(), V.rows());
	L.setZero();
	std::vector<Triplet<Scalar>> L_triplets;
	L_triplets.reserve(16 * F.rows());
	int simplex_size = F.cols();
	// 4 for tets
	assert(simplex_size == 4);

	// Loop over tetrahedra
	for (int i = 0; i < F.rows(); i++) {
		Matrix<Scalar, 3, 3> Dm, Dm_inv_trans;
		Dm << V(F(i, 0), 0) - V(F(i, 3), 0), V(F(i, 1), 0) - V(F(i, 3), 0), V(F(i, 2), 0) - V(F(i, 3), 0),
			  V(F(i, 0), 1) - V(F(i, 3), 1), V(F(i, 1), 1) - V(F(i, 3), 1), V(F(i, 2), 1) - V(F(i, 3), 1),
			  V(F(i, 0), 2) - V(F(i, 3), 2), V(F(i, 1), 2) - V(F(i, 3), 2), V(F(i, 2), 2) - V(F(i, 3), 2);
		Dm_inv_trans = Dm.inverse().transpose();
		
		//SparseMatrix<Scalar> Ai(3, V.rows());
		//vector<Triplet<Scalar>> Ai_IJV;
		//Ai_IJV.push_back(Triplet<Scalar>(0, F(i, 0), 1)); Ai_IJV.push_back(Triplet<Scalar>(0, F(i, 3), -1));
		//Ai_IJV.push_back(Triplet<Scalar>(1, F(i, 1), 1)); Ai_IJV.push_back(Triplet<Scalar>(1, F(i, 3), -1));
		//Ai_IJV.push_back(Triplet<Scalar>(2, F(i, 2), 1)); Ai_IJV.push_back(Triplet<Scalar>(2, F(i, 3), -1));
		//Ai.setFromTriplets(Ai_IJV.begin(), Ai_IJV.end());
		//SparseMatrix<Scalar> Gi(3, V.rows());
		//Gi = (Dm_inv_trans * Ai).sparseView();
		//L += Gi.transpose() * Gi;

		Eigen::Vector3d v4 = -Dm_inv_trans.col(0) - Dm_inv_trans.col(1) - Dm_inv_trans.col(2);
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {

				L_triplets.push_back(Triplet<Scalar>(
					F(i, row),
					F(i, col),
					Dm_inv_trans.col(row).dot(Dm_inv_trans.col(col))
					));
			}
			L_triplets.push_back(Triplet<Scalar>(
				F(i, row),
				F(i, 3),
				Dm_inv_trans.col(row).dot(v4)
				));
		}
		for (int col = 0; col < 3; col++) {
			L_triplets.push_back(Triplet<Scalar>(
					F(i, 3),
					F(i, col),
					Dm_inv_trans.col(col).dot(v4)
			));
		}
		L_triplets.push_back(Triplet<Scalar>(
			F(i, 3),
			F(i, 3),
			v4.dot(v4)
			));
	}
	L.setFromTriplets(L_triplets.begin(), L_triplets.end());
	L.makeCompressed();
}


template <typename DerivedV, typename DerivedF, typename Scalar>
IGL_INLINE void igl::laplacian_matrix(
	const Eigen::MatrixBase<DerivedV> & V,
	const Eigen::MatrixBase<DerivedF> & F,
	std::vector<Eigen::Triplet<Scalar>> & L_triplets)
{
	using namespace Eigen;
	using namespace std;

	L_triplets.reserve(16 * F.rows());
	int simplex_size = F.cols();
	// 4 for tets
	assert(simplex_size == 4);

	// Loop over tetrahedra
	for (int i = 0; i < F.rows(); i++) {
		Matrix<Scalar, 3, 3> Dm, Dm_inv_trans;
		Dm << V(F(i, 0), 0) - V(F(i, 3), 0), V(F(i, 1), 0) - V(F(i, 3), 0), V(F(i, 2), 0) - V(F(i, 3), 0),
			V(F(i, 0), 1) - V(F(i, 3), 1), V(F(i, 1), 1) - V(F(i, 3), 1), V(F(i, 2), 1) - V(F(i, 3), 1),
			V(F(i, 0), 2) - V(F(i, 3), 2), V(F(i, 1), 2) - V(F(i, 3), 2), V(F(i, 2), 2) - V(F(i, 3), 2);
		Dm_inv_trans = Dm.inverse().transpose();

		Eigen::Vector3d v4 = -Dm_inv_trans.col(0) - Dm_inv_trans.col(1) - Dm_inv_trans.col(2);
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {

				L_triplets.push_back(Triplet<Scalar>(
					F(i, row),
					F(i, col),
					Dm_inv_trans.col(row).dot(Dm_inv_trans.col(col))
					));
			}
			L_triplets.push_back(Triplet<Scalar>(
				F(i, row),
				F(i, 3),
				Dm_inv_trans.col(row).dot(v4)
				));
		}
		for (int col = 0; col < 3; col++) {
			L_triplets.push_back(Triplet<Scalar>(
				F(i, 3),
				F(i, col),
				Dm_inv_trans.col(col).dot(v4)
				));
		}
		L_triplets.push_back(Triplet<Scalar>(
			F(i, 3),
			F(i, 3),
			v4.dot(v4)
			));
	}
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::laplacian_matrix<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, Eigen::SparseMatrix<double, 0, int> &);
template void igl::laplacian_matrix< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::MatrixBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::vector< Eigen::Triplet<double, int>, std::allocator< Eigen::Triplet<double, int> > > &);

#endif
