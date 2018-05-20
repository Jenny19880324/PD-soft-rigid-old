// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "fit_rigid_motion.h"
#include "polar_svd3x3.h"
#include "polar_svd.h"
#include <vector>
#include <Eigen/Sparse>

template <typename DerivedV>
IGL_INLINE void igl::fit_rigid_motion(
	const Eigen::PlainObjectBase<DerivedV> & p,
	const Eigen::PlainObjectBase<DerivedV> & q,
	const Eigen::VectorXd & w,
	Eigen::Matrix3d & R,
	Eigen::RowVector3d & t)
{
	using namespace std;
	using namespace Eigen;
	int n = w.rows();
	
	// 1. Compute the weighted centroid of both point sets
	RowVector3d p_bar = RowVector3d::Zero();
	RowVector3d q_bar = RowVector3d::Zero();
	for (int i = 0; i < n; ++i) {
		p_bar += w(i) * p.row(i);
		q_bar += w(i) * q.row(i);
	}
	p_bar /= w.sum();
	q_bar /= w.sum();
	
	// 2. Compute the centered vectors
	MatrixXd X(n, 3);
	MatrixXd Y(n, 3);
	for (int i = 0; i < n; ++i) {
		X.row(i) = p.row(i) - p_bar;
		Y.row(i) = q.row(i) - q_bar;
	}
	
	// 3. Compute the dxd covariance matrix
	// S = XWY'
	typedef typename DerivedV::Scalar Scalar;
	vector<Triplet<Scalar>> triplets;
	for (int i = 0; i < n; ++i) {
		triplets.push_back(Triplet<Scalar>(i, i, w(i)));
	}
	SparseMatrix<double> W(n, n);
	W.setFromTriplets(triplets.begin(), triplets.end());
	
	Matrix3d S = X.transpose() * W * Y;
	
	// 4. Compute the sigular value decomposition S = U * Sigma * V'
	Matrix3d ti, ui, vi;
	Vector3d _;
	igl::polar_svd(S, R, ti, ui, _, vi);
	
	// 5. Compute the optimal translation as
	t = q_bar - p_bar * R;
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::fit_rigid_motion<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::Matrix<double, -1, 1, 0, -1, 1> const &, Eigen::Matrix<double, 3, 3, 0, 3, 3> &, Eigen::Matrix<double, 1, 3, 1, 1, 3> &);
#endif