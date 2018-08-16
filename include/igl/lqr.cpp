#include "lqr.h"
#include <iostream>
#include <Eigen/Dense>

template <
	typename DerivedX, 
	typename Scalar>
IGL_INLINE void igl::lqr(
	int N,
	Scalar dt,
	Scalar rho,
	const DerivedX &x_init,
	const DerivedX &x_target,
	std::vector<DerivedX> &x_sol)
{
	using namespace std;
	using namespace Eigen;

	Eigen::MatrixXd A, B, C;
	
	A.resize(6, 6);
	A << 1.0, 0.0, 0.0, dt, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0, dt, 0.0,
		0.0, 0.0, 1.0, 0.0, 0.0, dt,
		0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

	B.resize(6, 3);
	Scalar dt2 = 0.5 * dt * dt;
	B << dt2, 0.0, 0.0,
		0.0, dt2, 0.0,
		0.0, 0.0, dt2,
		dt, 0.0, 0.0,
		0.0, dt, 0.0,
		0.0, 0.0, dt;

	C.resize(3, 6);
	C << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0, 0.0, 0.0;

	Eigen::MatrixXd Q, Qf, R;
	Qf = C.transpose() * C;

	R = rho * Eigen::MatrixXd::Identity(3, 3);

	// follow Summary of LQR solution via DP
	// 1. set P_N:= Qf
	std::vector<Eigen::MatrixXd> P(N + 1);
	P[N] = Qf;

	// 2. for t = N, ..., 1,
	//    P_{t-1} = Q + A^T P_t A - A^T P_t B(R + B^T P_t * B)^(-1) B^T P_t A
	Eigen::MatrixXd A_trans = A.transpose();
	Eigen::MatrixXd B_trans = B.transpose();
	for (int t = N; t >= 1; t--) {
		Q = Qf;
		P[t - 1] = Q + A_trans * P[t] * A - A_trans * P[t] * B * (R + B_trans * P[t] * B).inverse() * B_trans * P[t] * A;
	}

	// 3. for t = 0, ..., N - 1, define K_t:= -(R + B_trans * P_t * B)^ (-1) * B_trans * P_{t+1} * A
	x_sol.resize(N);
	DerivedX x_t = x_init;
	for (int t = 0; t <= N - 1; t++) {
		Eigen::MatrixXd K_t = -(R + B_trans * P[t + 1] * B).inverse() * B_trans * P[t + 1] * A;
		Eigen::Vector3d u_t;

		u_t = K_t * (x_t - x_target);
		x_t = A * x_t + B * u_t;

		x_sol[t] = x_t;
	}
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::lqr<Eigen::Matrix<double, 6, 1, 0, 6, 1>, double>(int, double, double, Eigen::Matrix<double, 6, 1, 0, 6, 1> const &, Eigen::Matrix<double, 6, 1, 0, 6, 1> const &, std::vector<Eigen::Matrix<double, 6, 1, 0, 6, 1>, std::allocator<Eigen::Matrix<double, 6, 1, 0, 6, 1> > > &);

#endif
