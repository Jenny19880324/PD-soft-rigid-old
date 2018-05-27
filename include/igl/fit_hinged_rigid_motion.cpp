// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "fit_hinged_rigid_motion.h"
#include "fit_rigid_motion.h"
#include "min_quad_with_fixed.h"
#include "polar_svd3x3.h"
#include "polar_svd.h"
#include <vector>
#include <Eigen/Sparse>
#include <omp.h>


template <typename DerivedV>
IGL_INLINE void igl::fit_hinged_rigid_motion(
	const Eigen::PlainObjectBase<DerivedV> & v1,
	const Eigen::PlainObjectBase<DerivedV> & d1,
	const Eigen::PlainObjectBase<DerivedV> & cd1,
	const Eigen::PlainObjectBase<DerivedV> & v2,
	const Eigen::PlainObjectBase<DerivedV> & d2,
	const Eigen::PlainObjectBase<DerivedV> &cd2,
	const Eigen::RowVector3d &p,
	Eigen::Matrix3d &R1,
	Eigen::RowVector3d &t1,
	Eigen::Matrix3d &R2,
	Eigen::RowVector3d &t2
)
{
	using namespace std;
	using namespace Eigen;

	int max_iter = 1;
	int iter = 0;
	while (iter < max_iter) {
		// step 1. ignore the p2p constraint, find optimal R1_0 t1_0, R2_0, t2_0
		VectorXd weight_1 = VectorXd::Ones(v1.rows());
		VectorXd weight_2 = VectorXd::Ones(v2.rows());

		fit_rigid_motion(v1, cd1, weight_1, R1, t1);
		fit_rigid_motion(v2, cd2, weight_2, R2, t2);

		RowVector3d p1 = p * R1 + t1;
		RowVector3d p2 = p * R2 + t2;
		double obj = (p1 - p2).squaredNorm();

		std::cout << "obj = " << obj << std::endl;

		if (obj < 1e-2) {
			//return;
		}

		// construct A
		Matrix3d Q1tQ1 = Matrix3d::Zero();
		Matrix3d Q1tR1 = Matrix3d::Zero();
		Matrix3d R1tQ1 = Matrix3d::Zero();
		Matrix3d R1tR1 = Matrix3d::Zero();
		Matrix3d R1_trans = R1.transpose();
		RowVector3d Q1te1 = RowVector3d::Zero();
		RowVector3d R1te1 = RowVector3d::Zero();
		double e1te1 = 0.;
		for (int i = 0; i < v1.rows(); i++)
		{
			double x = v1.row(i).x();
			double y = v1.row(i).y();
			double z = v1.row(i).z();

			Matrix3d v_cross;
			v_cross << 0., -z, y,
				z, 0., -x,
				-y, x, 0.;
			Eigen::Matrix3d Q1 = R1_trans * v_cross;
			Eigen::Matrix3d Q1_trans = Q1.transpose();

			Q1tQ1 += Q1_trans * Q1;
			Q1tR1 += Q1_trans * R1_trans;
			R1tQ1 += R1 * Q1;
			R1tR1 += R1 * R1_trans;

			RowVector3d e1 = v1.row(i) * R1 + t1 - d1.row(i);
			Q1te1 += e1 * Q1;
			R1te1 += e1 * R1_trans;
			e1te1 += e1 * e1.transpose();
		}

		Matrix3d Q2tQ2 = Matrix3d::Zero();
		Matrix3d Q2tR2 = Matrix3d::Zero();
		Matrix3d R2tQ2 = Matrix3d::Zero();
		Matrix3d R2tR2 = Matrix3d::Zero();
		Matrix3d R2_trans = R2.transpose();
		RowVector3d Q2te2 = RowVector3d::Zero();
		RowVector3d R2te2 = RowVector3d::Zero();
		double e2te2 = 0.;
		for (int i = 0; i < v2.rows(); i++)
		{
			double x = v2.row(i).x();
			double y = v2.row(i).y();
			double z = v2.row(i).z();

			Matrix3d v_cross;
			v_cross << 0., -z, y,
				z, 0., -x,
				-y, x, 0.;
			Eigen::Matrix3d Q2 = R2_trans * v_cross;
			Eigen::Matrix3d Q2_trans = Q2.transpose();

			Q2tQ2 += Q2_trans * Q2;
			Q2tR2 += Q2_trans * R2_trans;
			R2tQ2 += R2 * Q2;
			R2tR2 += R2 * R2_trans;

			RowVector3d e2 = v2.row(i) * R2 + t2 - d2.row(i);
			Q2te2 += e2 * Q2;
			R2te2 += e2 * R2_trans;
			e2te2 += e2 * e2.transpose();
		}

		MatrixXd A;
		A.resize(12, 12);

		A.setZero();
		A.block(0, 0, 3, 3) = Q1tQ1;
		A.block(0, 3, 3, 3) = -Q1tR1;
		A.block(3, 0, 3, 3) = -R1tQ1;
		A.block(3, 3, 3, 3) = R1tR1;

		A.block(6, 6, 3, 3) = Q2tQ2;
		A.block(6, 9, 3, 3) = -Q2tR2;
		A.block(9, 6, 3, 3) = -R2tQ2;
		A.block(9, 9, 3, 3) = R2tR2;
		SparseMatrix<double> A_s = A.sparseView();

		// construct B
		VectorXd B;
		B.resize(12, 1);
		B.setZero();
		B.block(0, 0, 3, 1) = -Q1te1.transpose();
		B.block(3, 0, 3, 1) = R1te1.transpose();
		B.block(6, 0, 3, 1) = -Q2te2.transpose();
		B.block(9, 0, 3, 1) = R2te2.transpose();

		// construct Aeq
		MatrixXd Aeq;
		Aeq.resize(3, 12);
		Aeq.setZero();
		Matrix3d p_cross;
		p_cross << 0., -p.z(), p.y(),
			p.z(), 0., -p.x(),
			-p.y(), p.x(), 0.;
		Aeq.block(0, 0, 3, 3) = R1_trans * p_cross;
		Aeq.block(0, 3, 3, 3) = -R1_trans;
		Aeq.block(0, 6, 3, 3) = -R2_trans * p_cross;
		Aeq.block(0, 9, 3, 3) = R2_trans;
		Eigen::SparseMatrix<double> Aeq_s = Aeq.sparseView();

		// construct Beq
		MatrixXd Beq;
		Beq.resize(3, 1);
		Beq.setZero();
		Beq = (p * (R1 - R2) + t1 - t2).transpose();

		VectorXd Z;
		min_quad_with_fixed_data<double> solver_data;
		min_quad_with_fixed_precompute(A_s, VectorXi(), Aeq_s, true, solver_data);
		min_quad_with_fixed_solve(solver_data, B, MatrixX3d(), Beq, Z);

		Eigen::Vector3d w1 = Z.block(0, 0, 3, 1);
		Eigen::Vector3d l1 = Z.block(3, 0, 3, 1);
		Eigen::Vector3d w2 = Z.block(6, 0, 3, 1);
		Eigen::Vector3d l2 = Z.block(9, 0, 3, 1);

		Eigen::Matrix3d w1_cross;
		Eigen::Matrix3d w2_cross;
		w1_cross << 0., -w1.z(), w1.y(),
			w1.z(), 0., -w1.x(),
			-w1.y(), w1.x(), 0.;

		w2_cross << 0., -w2.z(), w2.y(),
			w2.z(), 0., -w2.x(),
			-w2.y(), w2.x(), 0.;

		std::cout << "t1 = " << t1 << std::endl;
		R1 = (R1_trans * (Matrix3d::Identity() + w1_cross)).transpose();
		t1 = (R1_trans * l1).transpose() + t1;
		std::cout << "t1 = " << t1 << std::endl;

		std::cout << "t2 = " << t2 << std::endl;
		R2 = (R2_trans * (Matrix3d::Identity() + w2_cross)).transpose();
		t2 = (R2_trans * l2).transpose() + t2;
		std::cout << "t2 = " << t2 << std::endl;

		Eigen::MatrixX3d td1 = v1 * R1 + t1.replicate(v1.rows(), 1);
		Eigen::MatrixX3d td2 = v2 * R2 + t2.replicate(v2.rows(), 1);
		double dist_1 = (td1 - d1).squaredNorm() + (td2 - d2).squaredNorm();
		std::cout << "dist_1 = " << dist_1 << std::endl;

		iter++;
	}

}

	template <typename DerivedV>
	IGL_INLINE void igl::fit_hinged_rigid_motion(
		const Eigen::PlainObjectBase<DerivedV> & v1,
		const Eigen::PlainObjectBase<DerivedV> & d1,
		const Eigen::PlainObjectBase<DerivedV> & v2,
		const Eigen::PlainObjectBase<DerivedV> & d2,
		const Eigen::RowVector3d &p,
		Eigen::Matrix3d &R1,
		Eigen::RowVector3d &t1,
		Eigen::Matrix3d &R2,
		Eigen::RowVector3d &t2
	)
	{
		using namespace std;
		using namespace Eigen;

		int max_iter = 10;
		int iter = 0;
		Eigen::MatrixXd cd1 = d1;
		Eigen::MatrixXd cd2 = d2;
		while (iter < max_iter) {
			// step 1. ignore the p2p constraint, find optimal R1_0 t1_0, R2_0, t2_0
			VectorXd weight_1 = VectorXd::Ones(v1.rows());
			VectorXd weight_2 = VectorXd::Ones(v2.rows());

			fit_rigid_motion(v1, cd1, weight_1, R1, t1);
			fit_rigid_motion(v2, cd2, weight_2, R2, t2);

			// construct A
			Matrix3d Q1tQ1 = Matrix3d::Zero();
			Matrix3d Q1tR1 = Matrix3d::Zero();
			Matrix3d R1tQ1 = Matrix3d::Zero();
			Matrix3d R1tR1 = Matrix3d::Zero();
			Matrix3d R1_trans = R1.transpose();
			RowVector3d Q1te1 = RowVector3d::Zero();
			RowVector3d R1te1 = RowVector3d::Zero();
			double e1te1 = 0.;
			for (int i = 0; i < v1.rows(); i++)
			{
				double x = v1.row(i).x();
				double y = v1.row(i).y();
				double z = v1.row(i).z();

				Matrix3d v_cross;
				v_cross << 0., -z, y,
					z, 0., -x,
					-y, x, 0.;
				Eigen::Matrix3d Q1 = R1_trans * v_cross;
				Eigen::Matrix3d Q1_trans = Q1.transpose();

				Q1tQ1 += Q1_trans * Q1;
				Q1tR1 += Q1_trans * R1_trans;
				R1tQ1 += R1 * Q1;
				R1tR1 += R1 * R1_trans;

				RowVector3d e1 = v1.row(i) * R1 + t1 - d1.row(i);
				Q1te1 += e1 * Q1;
				R1te1 += e1 * R1_trans;
				e1te1 += e1 * e1.transpose();
			}

			Matrix3d Q2tQ2 = Matrix3d::Zero();
			Matrix3d Q2tR2 = Matrix3d::Zero();
			Matrix3d R2tQ2 = Matrix3d::Zero();
			Matrix3d R2tR2 = Matrix3d::Zero();
			Matrix3d R2_trans = R2.transpose();
			RowVector3d Q2te2 = RowVector3d::Zero();
			RowVector3d R2te2 = RowVector3d::Zero();
			double e2te2 = 0.;
			for (int i = 0; i < v2.rows(); i++)
			{
				double x = v2.row(i).x();
				double y = v2.row(i).y();
				double z = v2.row(i).z();

				Matrix3d v_cross;
				v_cross << 0., -z, y,
					z, 0., -x,
					-y, x, 0.;
				Eigen::Matrix3d Q2 = R2_trans * v_cross;
				Eigen::Matrix3d Q2_trans = Q2.transpose();

				Q2tQ2 += Q2_trans * Q2;
				Q2tR2 += Q2_trans * R2_trans;
				R2tQ2 += R2 * Q2;
				R2tR2 += R2 * R2_trans;

				RowVector3d e2 = v2.row(i) * R2 + t2 - d2.row(i);
				Q2te2 += e2 * Q2;
				R2te2 += e2 * R2_trans;
				e2te2 += e2 * e2.transpose();
			}

			MatrixXd A;
			A.resize(12, 12);

			A.setZero();
			A.block(0, 0, 3, 3) = Q1tQ1;
			A.block(0, 3, 3, 3) = -Q1tR1;
			A.block(3, 0, 3, 3) = -R1tQ1;
			A.block(3, 3, 3, 3) = R1tR1;

			A.block(6, 6, 3, 3) = Q2tQ2;
			A.block(6, 9, 3, 3) = -Q2tR2;
			A.block(9, 6, 3, 3) = -R2tQ2;
			A.block(9, 9, 3, 3) = R2tR2;
			SparseMatrix<double> A_s = A.sparseView();

			// construct B
			VectorXd B;
			B.resize(12, 1);
			B.setZero();
			B.block(0, 0, 3, 1) = -Q1te1.transpose();
			B.block(3, 0, 3, 1) = R1te1.transpose();
			B.block(6, 0, 3, 1) = -Q2te2.transpose();
			B.block(9, 0, 3, 1) = R2te2.transpose();

			// construct Aeq
			MatrixXd Aeq;
			Aeq.resize(3, 12);
			Aeq.setZero();
			Matrix3d p_cross;
			p_cross << 0., -p.z(), p.y(),
				p.z(), 0., -p.x(),
				-p.y(), p.x(), 0.;
			Aeq.block(0, 0, 3, 3) = R1_trans * p_cross;
			Aeq.block(0, 3, 3, 3) = -R1_trans;
			Aeq.block(0, 6, 3, 3) = -R2_trans * p_cross;
			Aeq.block(0, 9, 3, 3) = R2_trans;
			Eigen::SparseMatrix<double> Aeq_s = Aeq.sparseView();

			// construct Beq
			MatrixXd Beq;
			Beq.resize(3, 1);
			Beq.setZero();
			Beq = (p * (R1 - R2) + t1 - t2).transpose();

			VectorXd Z;
			min_quad_with_fixed_data<double> solver_data;
			min_quad_with_fixed_precompute(A_s, VectorXi(), Aeq_s, true, solver_data);
			min_quad_with_fixed_solve(solver_data, B, VectorXd(), Beq, Z);

			Eigen::Vector3d w1 = Z.block(0, 0, 3, 1);
			Eigen::Vector3d l1 = Z.block(3, 0, 3, 1);
			Eigen::Vector3d w2 = Z.block(6, 0, 3, 1);
			Eigen::Vector3d l2 = Z.block(9, 0, 3, 1);

			Eigen::Matrix3d w1_cross;
			Eigen::Matrix3d w2_cross;
			w1_cross << 0., -w1.z(), w1.y(),
				w1.z(), 0., -w1.x(),
				-w1.y(), w1.x(), 0.;

			w2_cross << 0., -w2.z(), w2.y(),
				w2.z(), 0., -w2.x(),
				-w2.y(), w2.x(), 0.;

			R1 = (R1_trans * (Matrix3d::Identity() + w1_cross)).transpose();
			t1 = (R1_trans * l1).transpose() + t1;

			R2 = (R2_trans * (Matrix3d::Identity() + w2_cross)).transpose();
			t2 = (R2_trans * l2).transpose() + t2;

			Eigen::MatrixX3d cd1 = v1 * R1 + t1.replicate(v1.rows(), 1);
			Eigen::MatrixX3d cd2 = v2 * R2 + t2.replicate(v2.rows(), 1);

			RowVector3d p1 = p * R1 + t1;
			RowVector3d p2 = p * R2 + t2;
			double obj = (p1 - p2).squaredNorm();

			double dist_1 = (cd1 - d1).squaredNorm() + (cd2 - d2).squaredNorm();
			std::cout << "dist_1 = " << dist_1 << std::endl;

			if (obj < 1e-2) {
				return;
			}
			iter++;
		}
}

template <typename DerivedV, typename DerivedN>
IGL_INLINE void igl::fit_hinged_rigid_motion(
	const Eigen::PlainObjectBase<DerivedV> & V,
	const Eigen::PlainObjectBase<DerivedN> & N,
	const Eigen::PlainObjectBase<DerivedV> & P,
	const std::vector<std::vector<int>> & I,
	Eigen::PlainObjectBase<DerivedV> & U)
{
	using namespace std;
	using namespace Eigen;

	int max_iter = 10;
	int iter = 0;
	const int dim = 3;

	int number_of_rigid_bodies = N.rows() - 1;
	int number_of_joints = P.rows();
	int number_of_constraints = 0;
	for (int i = 0; i < I.size(); i++) {
		number_of_constraints += I[i].size() - 1;
	}
	//std::cout << "number_of_constraints = " << number_of_constraints << std::endl;
	VectorXi acum_N(number_of_rigid_bodies);
	acum_N(0) = N(0);
	for (int b_i = 1; b_i < number_of_rigid_bodies; b_i++) {
		acum_N(b_i) = acum_N(b_i - 1) + N(b_i);
	}
	assert(P.rows() == I.size());
	MatrixXd cU = U;

	MatrixXd A, Aeq;
	VectorXd B, Beq;
	A.resize(6 * number_of_rigid_bodies, 6 * number_of_rigid_bodies);
	Aeq.resize(3 * number_of_constraints, 6 * number_of_rigid_bodies);
	B.resize(6 * number_of_rigid_bodies, 1);
	Beq.resize(3 * number_of_constraints, 1);
	A.setZero(); Aeq.setZero();
	B.setZero(); Beq.setZero();

	MatrixX3d R, R_trans; // number_of_rigid_bodies * 3 by dim
	R.resize(3 * number_of_rigid_bodies, Eigen::NoChange);
	R_trans.resize(3 * number_of_rigid_bodies, Eigen::NoChange);
	MatrixX3d t; // number_of_rigid_bodies by dim
	t.resize(number_of_rigid_bodies, Eigen::NoChange);

	while (iter < max_iter) {
		#pragma omp parallel for
		for (int b_i = 0; b_i < number_of_rigid_bodies; b_i++) {
			MatrixXd cUb = cU.block(acum_N(b_i), 0, N(b_i + 1), dim);
			MatrixXd Ub = U.block(acum_N(b_i), 0, N(b_i + 1), dim);
			MatrixXd Vb = V.block(acum_N(b_i), 0, N(b_i + 1), dim);
			VectorXd weight = VectorXd::Ones(N(b_i + 1));
			Matrix3d Rb;
			RowVector3d tb;
			fit_rigid_motion(Vb, cUb, weight, Rb, tb);
			
			Matrix3d QbtQb = Matrix3d::Zero();
			Matrix3d QbtRb = Matrix3d::Zero();
			Matrix3d RbtQb = Matrix3d::Zero();
			Matrix3d RbtRb = Matrix3d::Zero();
			Matrix3d Rb_trans = Rb.transpose();
			RowVector3d Qbteb = RowVector3d::Zero();
			RowVector3d Rbteb = RowVector3d::Zero();
			double ebteb = 0.;
			R_trans.block(3 * b_i, 0, 3, 3) = Rb_trans;
			R.block(3 * b_i, 0, 3, 3) = Rb;
			t.block(b_i, 0, 1, 3) = tb;

			for (int i = 0; i < N(b_i + 1); i++) {
				double x = Vb.row(i).x();
				double y = Vb.row(i).y();
				double z = Vb.row(i).z();

				Matrix3d v_cross;
				v_cross << 0., -z, y,
					z, 0., -x,
					-y, x, 0.;

				Matrix3d Qb = Rb_trans * v_cross;
				Matrix3d Qb_trans = Qb.transpose();

				QbtQb += Qb_trans * Qb;
				QbtRb += Qb_trans * Rb_trans;
				RbtQb += Rb * Qb;
				RbtRb += Rb * Rb_trans;

				RowVector3d eb = Vb.row(i) * Rb + tb - Ub.row(i);
				Qbteb += eb * Qb;
				Rbteb += eb * Rb_trans;
				ebteb += eb * eb.transpose();

				// construct A 
				// todo use triplets, might be faster
				A.block(6 * b_i    , 6 * b_i    , 3, 3) = QbtQb;
				A.block(6 * b_i    , 6 * b_i + 3, 3, 3) = -QbtRb;
				A.block(6 * b_i + 3, 6 * b_i    , 3, 3) = -RbtQb;
				A.block(6 * b_i + 3, 6 * b_i + 3, 3, 3) = RbtRb;

				// construct B
				B.block(6 * b_i    , 0, 3, 1) = -Qbteb.transpose();
				B.block(6 * b_i + 3, 0, 3, 1) = Rbteb.transpose();
			}
		}
		
		#pragma omp parallel for
		for (int j_i = 0; j_i < number_of_joints; j_i++)
		{
			const RowVector3d p = P.row(j_i);
			double x = P.row(j_i).x();
			double y = P.row(j_i).y();
			double z = P.row(j_i).z();

			Matrix3d p_cross;
			p_cross << 0., -z, y,
						z, 0., -x,
						-y, x, 0.;

			int idx_1 = I[j_i][0];
			Matrix3d R1_trans = R_trans.block(idx_1 * 3, 0, 3, 3);
			Matrix3d R1 = R.block(idx_1 * 3, 0, 3, 3);
			RowVector3d t1 = t.block(idx_1, 0, 1, 3);
			int constraint_idx = 0;
			for (int i = 0; i < j_i; i++) {
				constraint_idx += I[j_i].size() - 1;
			}
			Aeq.block(3 * constraint_idx, 6 * idx_1    , 3, 3) = R1_trans * p_cross;
			Aeq.block(3 * constraint_idx, 6 * idx_1 + 3, 3, 3) = -R1_trans;
			for (int i = 1; i < I[j_i].size(); i++) {
				int idx_2 = I[j_i][i];
				Matrix3d R2_trans = R_trans.block(idx_2 * 3, 0, 3, 3);
				Matrix3d R2 = R.block(idx_2 * 3, 0, 3, 3);
				RowVector3d t2 = t.block(idx_2, 0, 1, 3);
				Aeq.block(3 * constraint_idx, 6 * idx_2    , 3, 3) = -R2_trans * p_cross;
				Aeq.block(3 * constraint_idx, 6 * idx_2 + 3, 3, 3) = R2_trans;
				Beq.block(3 * constraint_idx, 0, 3, 1) = (p * (R1 - R2) + t1 - t2).transpose();
				constraint_idx++;
			}
		}

		SparseMatrix<double> A_s = A.sparseView();
		SparseMatrix<double> Aeq_s = Aeq.sparseView();

		VectorXd Z;
		min_quad_with_fixed_data<double> solver_data;
		min_quad_with_fixed_precompute(A_s, VectorXi(), Aeq_s, true, solver_data);
		min_quad_with_fixed_solve(solver_data, B, VectorXd(), Beq, Z);

		double dist = 0.;
		#pragma omp parallel for
		for (int b_i = 0; b_i < number_of_rigid_bodies; b_i++) {
			Vector3d wb = Z.block(b_i * 6    , 0, 3, 1);
			Vector3d lb = Z.block(b_i * 6 + 3, 0, 3, 1);

			Matrix3d wb_cross;
			wb_cross << 0., -wb.z(), wb.y(),
				wb.z(), 0., -wb.x(),
				-wb.y(), wb.x(), 0.;

			Matrix3d Rb_trans = R_trans.block(b_i * 3, 0, 3, 3);
			Matrix3d Rb = (Rb_trans * (Matrix3d::Identity() + wb_cross)).transpose();
			RowVector3d tb = (Rb_trans * lb).transpose() + t.block(b_i, 0, 1, 3);
			R.block(b_i * 3, 0, 3, 3) = Rb;
			R_trans.block(b_i * 3, 0, 3, 3) = Rb.transpose();

			MatrixXd Vb = V.block(acum_N(b_i), 0, N(b_i + 1), dim);
			MatrixXd cUb = Vb * Rb + tb.replicate(N(b_i + 1), 1);
			MatrixXd Ub = U.block(acum_N(b_i), 0, N(b_i + 1), dim);
		    cU.block(acum_N(b_i), 0, N(b_i + 1), dim) = cUb;
			dist += (cUb - Ub).squaredNorm();
		}
		//std::cout << "dist = " << dist << std::endl;

		double obj = 0.;
		for (int j_i = 0; j_i < number_of_joints; j_i++)
		{
			const RowVector3d p = P.row(j_i);
			int idx_1 = I[j_i][0];
			Matrix3d R1 = R.block(idx_1 * 3, 0, 3, 3);
			RowVector3d t1 = t.block(idx_1, 0, 1, 3);
			RowVector3d p1 = p * R1 + t1;
			for (int i = 1; i < I[j_i].size(); i++) {
				int idx_2 = I[j_i][i];
				Matrix3d R2 = R.block(idx_2 * 3, 0, 3, 3);
				RowVector3d t2 = t.block(idx_2, 0, 1, 3);
				RowVector3d  p2 = p * R2 + t2;
				obj += (p1 - p2).squaredNorm();
			}
		}
		//std::cout << "obj = " << obj << std::endl;

		if (obj < 1e-6) {
			break;
		}
		iter++;
	}
	U = cU;
}


#ifdef IGL_STATIC_LIBRARY
template void igl::fit_hinged_rigid_motion<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::Matrix<double, 1, 3, 1, 1, 3> const &, Eigen::Matrix<double, 3, 3, 0, 3, 3> &, Eigen::Matrix<double, 1, 3, 1, 1, 3> &, Eigen::Matrix<double, 3, 3, 0, 3, 3> &, Eigen::Matrix<double, 1, 3, 1, 1, 3> &);
template void igl::fit_hinged_rigid_motion<Eigen::Matrix<double,-1,-1,0,-1,-1> >(Eigen::PlainObjectBase<Eigen::Matrix<double,-1,-1,0,-1,-1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double,-1,-1,0,-1,-1> > const &,Eigen::PlainObjectBase<Eigen::Matrix<double,-1,-1,0,-1,-1> > const &,Eigen::PlainObjectBase<Eigen::Matrix<double,-1,-1,0,-1,-1> > const &,Eigen::Matrix<double,1,3,1,1,3> const &,Eigen::Matrix<double,3,3,0,3,3> &,Eigen::Matrix<double,1,3,1,1,3> &,Eigen::Matrix<double,3,3,0,3,3> &,Eigen::Matrix<double,1,3,1,1,3> &);
template void igl::fit_hinged_rigid_motion<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, std::vector<std::vector<int, std::allocator<int> >, std::allocator< std::vector<int, std::allocator<int> > > > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);
//template void igl::fit_hinged_rigid_motion<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);
#endif