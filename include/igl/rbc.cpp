#include "rbc.h"
#include "shape_matrix.h"
#include "massmatrix.h"
#include "fit_rotations.h"
#include "fit_rigid_motion.h"
#include "fit_hinged_rigid_motion.h"
#include "laplacian_matrix.h"
#include "rbc_rhs.h"
#include "self_collision.h"


extern Eigen::MatrixXd P;
extern Eigen::MatrixXi SF;
extern Eigen::VectorXi SV;
extern std::vector<std::vector<int>> I;

template<
	typename DerivedV,
	typename DerivedF,
	typename DerivedN,
	typename Derivedb>
IGL_INLINE bool igl::rbc_precomputation(
	const Eigen::PlainObjectBase<DerivedV> & V,
	const Eigen::PlainObjectBase<DerivedF> & F,
	const Eigen::PlainObjectBase<DerivedN> & N,
	const int dim,
	const Eigen::PlainObjectBase<Derivedb> & b,
	RBCData & data)
{
	using namespace std;
	using namespace Eigen;
	typedef typename DerivedV::Scalar Scalar;
	// number of vertices
	const int n = V.rows();

	data.n = n;
	data.N = N;
	assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds");
	assert((b.size() == 0 || b.minCoeff() >=0) && "b out of bounds");
	// remember b
	data.b = b;
	data.F = F;
	data.V = V;
	// dimension
	assert((dim == 3 || dim == 2) && "dim should be 2 or 3");
	data.dim = dim;
	// Defaults
	//data.f_ext = MatrixXd::Zero(n, data.dim);
	
	assert(data.dim <= V.cols() && "solve dim should be <= embedding");
	
	RBCEnergyType eff_energy = data.energy;
	if(eff_energy == RBC_ENERGY_TYPE_DEFAULT)
	{
		switch(F.cols())
		{
			case 3:
				break;
			case 4:
				eff_energy = RBC_ENERGY_TYPE_RBC;
				break;
			default:
				assert(false);
		}
	}
	
	shape_matrix(V, F, eff_energy, data.SM);
	typedef SparseMatrix<Scalar> SparseMatrixS;
	SparseMatrixS L;

	laplacian_matrix(V, F, L);
	
	rbc_rhs(V, F, data.dim, eff_energy, data.J);

	if (data.constraint == SOFT_CONSTRAINT) {
		for (int i = 0; i < b.rows(); i++) {
			SparseMatrix<Scalar> Ai(1, V.rows());
			vector<Triplet<Scalar>> Ai_IJV;
			Ai_IJV.push_back(Triplet<Scalar>(0, b(i), 1.0));
			Ai.setFromTriplets(Ai_IJV.begin(), Ai_IJV.end());
			L += data.constraint_weight * Ai.transpose() * Ai;
		}
	}

	SparseMatrix<double> Q = L.eval();

	if (data.with_dynamics)
	{
		const double h = data.h;
		assert(h != 0);
		SparseMatrix<double> M;
		massmatrix(V, F, MASSMATRIX_TYPE_DEFAULT, data.M);
		data.M *= data.mass_scaling;
		const double dw = 1. / data.mu;
		SparseMatrix<double> DQ = dw * 1./(h * h) * data.M;
		Q += DQ;
		// Dummy external forces
		//data.f_ext = MatrixXd::Zero(n, data.dim);
		data.f_ext = Eigen::RowVector3d(0.,(double)data.g, 0.).replicate(V.rows(), 1);
		data.f_ext = data.M * data.f_ext;
		if (data.vel.rows() == 0) {
			data.vel = MatrixXd::Zero(n, data.dim);
		}
	}

	if (eff_energy == RBC_ENERGY_TYPE_RBC) {
		data.nf = N(0);
		data.nb = 0;
		for (int i = 1; i < N.rows(); i++) {
			data.nb += N(i);
		}
		assert(data.n == data.nf + data.nb);

		// construct Ab
		data.m = N.rows() - 1; // N(0) is the number of vertices of the elastic part
		assert(data.m > 0);
		data.Ab.resize(data.nb, 4 * data.m);
		data.Ab.setZero();
		int row = 0;
		for (int i = 1; i < N.rows(); i++)
		{
			data.Ab.block(row, (i - 1) * 4, N(i), 4) << V.block(data.nf + row, 0, N(i), 3), MatrixXd::Constant(N(i), 1, 1);
			row += N(i);
		}
		assert(row == data.nb);

		data.B.resize(data.n, data.nf + 4 * data.m);
		std::vector<Triplet<double>> B_IJV;
		for (int i = 0; i < data.nf; i++) {
			B_IJV.push_back(Triplet<double>(i, i, 1));
		}
		for (int i = 0; i < data.nb; i++) {
			for (int j = 0; j < 4 * data.m; j++) {
				B_IJV.push_back(Triplet<double>(data.nf + i, data.nf + j, data.Ab(i, j)));
			}
		}
		data.B.setFromTriplets(B_IJV.begin(), B_IJV.end());
		data.B_trans = data.B.transpose();
		Q = data.B_trans * Q * data.B;
	}

	if (data.constraint == SOFT_CONSTRAINT) {
		return min_quad_with_fixed_precompute(
			Q, VectorXi(), SparseMatrix<double>(), true, data.solver_data);
	}

	return min_quad_with_fixed_precompute(
		Q, b, SparseMatrix<double>(), true, data.solver_data);
}


// when collision added, A is altered each step
template <
	typename DerivedV,
	typename DerivedF,
	typename DerivedN,
	typename Derivedb>
	IGL_INLINE bool igl::rbc_precomputation(
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedV> & U,
		const Eigen::PlainObjectBase<DerivedF> & F,
		const Eigen::PlainObjectBase<DerivedN> & N,
		const int dim,
		const Eigen::PlainObjectBase<Derivedb> & b,
		RBCData & data)
{

	using namespace std;
	using namespace Eigen;
	typedef typename DerivedV::Scalar Scalar;
	// number of vertices
	const int n = V.rows();

	data.n = n;
	data.N = N;
	assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds");
	assert((b.size() == 0 || b.minCoeff() >= 0) && "b out of bounds");
	// remember b
	data.b = b;
	data.F = F;
	data.V = V;
	// dimension
	assert((dim == 3 || dim == 2) && "dim should be 2 or 3");
	data.dim = dim;
	// Defaults
	//data.f_ext = MatrixXd::Zero(n, data.dim);

	assert(data.dim <= V.cols() && "solve dim should be <= embedding");

	RBCEnergyType eff_energy = data.energy;
	if (eff_energy == RBC_ENERGY_TYPE_DEFAULT)
	{
		switch (F.cols())
		{
		case 3:
			break;
		case 4:
			eff_energy = RBC_ENERGY_TYPE_RBC;
			break;
		default:
			assert(false);
		}
	}

	shape_matrix(V, F, eff_energy, data.SM);
	typedef SparseMatrix<Scalar> SparseMatrixS;
	SparseMatrixS L;
	std::vector<Triplet<Scalar>> L_triplets;

	//laplacian_matrix(V, F, L);
	laplacian_matrix(V, F, L_triplets);

	if (data.constraint == SOFT_CONSTRAINT) {
		for (int i = 0; i < b.rows(); i++) {
			L_triplets.push_back(Triplet<Scalar>(b(i), b(i), data.constraint_weight));
		}
	}


	if (data.collision_enabled) {
		double floor_y = data.floor_y;
		for (int i = 0; i < U.rows(); i++) {
			if (U.row(i).y() < floor_y) {
				L_triplets.push_back(Triplet<Scalar>(i, i, data.collision_weight));
			}
		}
	}

	L.resize(V.rows(), V.rows());
	L.setZero();
	L.setFromTriplets(L_triplets.begin(), L_triplets.end());
	L.makeCompressed();
	SparseMatrix<double> Q = L.eval();

	if (data.with_dynamics)
	{
		const double h = data.h;
		assert(h != 0);
		SparseMatrix<double> M;
		massmatrix(V, F, MASSMATRIX_TYPE_DEFAULT, data.M);
		data.M *= data.mass_scaling;
		const double dw = 1. / data.mu;
		SparseMatrix<double> DQ = dw * 1. / (h * h) * data.M;
		Q += DQ;
		// Dummy external forces
		//data.f_ext = MatrixXd::Zero(n, data.dim);
		data.f_ext = Eigen::RowVector3d(0., (double)data.g, 0.).replicate(V.rows(), 1);
		data.f_ext = data.M * data.f_ext;	
	}

	if (eff_energy == RBC_ENERGY_TYPE_RBC) {
		data.nf = N(0);
		data.nb = 0;
		for (int i = 1; i < N.rows(); i++) {
			data.nb += N(i);
		}
		assert(data.n == data.nf + data.nb);

		// construct Ab
		data.m = N.rows() - 1; // N(0) is the number of vertices of the elastic part
		assert(data.m > 0);
		data.Ab.resize(data.nb, 4 * data.m);
		data.Ab.setZero();
		int row = 0;
		for (int i = 1; i < N.rows(); i++)
		{
			data.Ab.block(row, (i - 1) * 4, N(i), 4) << V.block(data.nf + row, 0, N(i), 3), MatrixXd::Constant(N(i), 1, 1);
			row += N(i);
		}
		assert(row == data.nb);

		data.B.resize(data.n, data.nf + 4 * data.m);
		std::vector<Triplet<double>> B_IJV;
		for (int i = 0; i < data.nf; i++) {
			B_IJV.push_back(Triplet<double>(i, i, 1));
		}
		for (int i = 0; i < data.nb; i++) {
			for (int j = 0; j < 4 * data.m; j++) {
				B_IJV.push_back(Triplet<double>(data.nf + i, data.nf + j, data.Ab(i, j)));
			}
		}
		data.B.setFromTriplets(B_IJV.begin(), B_IJV.end());
		data.B_trans = data.B.transpose();
		Q = data.B_trans * Q * data.B;
	}

	if (data.constraint == SOFT_CONSTRAINT) {
		return min_quad_with_fixed_precompute(
			Q, VectorXi(), SparseMatrix<double>(), true, data.solver_data);
	}

	return min_quad_with_fixed_precompute(
		Q, b, SparseMatrix<double>(), true, data.solver_data);
}



template <
	typename Derivedbc,
	typename DerivedU>
IGL_INLINE bool igl::rbc_solve(
	const Eigen::PlainObjectBase<Derivedbc> &bc,
	RBCData &data,
	Eigen::PlainObjectBase<DerivedU> & U)
	{
		using namespace Eigen;
		using namespace std;
 		if (data.b.size() != bc.rows()) {
			std::cout << "data.b.size() = " << data.b.size() << std::endl;
			std::cout << "bc.rows()     = " << bc.rows() << std::endl;
		}
		assert(data.b.size() == bc.rows());
		if (bc.size() > 0)
		{
			assert(bc.cols() == data.dim && "bc.cols() match data.dim");
		}
		const int n = data.n;
		int iter = 0;
		if (U.size() == 0)
		{
			// terrible initial guess.. should at least copy input mesh
#ifndef NDEBUG
			cerr << "rbc_solve: Using terrible initial guess for U. Try U = V." << endl;
#endif
			U = MatrixXd::Zero(data.n, data.dim);
		}else
		{
			assert(U.cols() == data.dim && "U.cols() match data.dim");
		}
		// changes each rbc iteration
		MatrixXd U_prev = U;
		// doesn't change for fixed with_dynamics timestep
		MatrixXd U0;
		if (data.with_dynamics)
		{
			U0 = U_prev;
		}
		while(iter < data.max_iter)
		{
			U_prev = U;
			// enforce boundary conditions exactly
			for (int bi = 0; bi < bc.rows(); bi++)
			{
				U.row(data.b(bi)) = bc.row(bi);
			}
			
			const int Rdim = data.dim;
			MatrixXd R(data.J.cols(), Rdim);
			MatrixXd S(data.SM.rows(), Rdim);
			MatrixXd F(data.SM.rows(), Rdim);

			shape_matrix(U, data.F, data.energy, S);
			const int nr = data.SM.rows() / 3;
			for (int r = 0; r < nr; ++r) {
				Matrix3d Dm = data.SM.block(3 * r, 0, 3, 3).transpose();
				Matrix3d Ds = S.block(3 * r, 0, 3, 3).transpose();
				F.block(3 * r, 0, 3, 3) = Ds * Dm.inverse();
			}
			fit_rotations(F, R);

			MatrixXd Dl;
			if (data.with_dynamics)
			{ 
				assert(data.M.rows() == n && 
				"No mass matrix. Call rbc_precomputation if chaning with_dynamics");
				const double h = data.h;
				assert(h != 0);
				const double dw = 1. / data.mu;
				Dl = dw * (1./(h*h)*data.M*(-U0-h*data.vel) - data.f_ext);
			}
			
			Matrix<double, Eigen::Dynamic, 3> B = -data.J * R;
			Matrix<double, Eigen::Dynamic, 3> Beq;
			Matrix<double, Eigen::Dynamic, 3> q;
			assert(B.rows() == data.n);
			assert(B.cols() == data.dim);
			
			if (data.with_dynamics)
			{
				B += Dl;
			}

			if (data.collision_enabled) {
				double floor_y = data.floor_y;
				for (int i = 0; i < U0.rows(); i++) {
					if (U0.row(i).y() < floor_y) {
						double x = U0.row(i).x();
						double y = floor_y;
						double z = U0.row(i).z();
						B.row(i) -= data.collision_weight * Eigen::RowVector3d(x, floor_y, z);
					}
				}
			}

			if (data.constraint == SOFT_CONSTRAINT) {
				for (int i = 0; i < data.b.rows(); i++) {
					B.row(data.b(i)) -= data.constraint_weight * bc.row(i);
				}
			}



			if (data.energy == RBC_ENERGY_TYPE_RBC) {
				B = data.B_trans * B;
			}
			
			if (data.constraint == SOFT_CONSTRAINT) {
				min_quad_with_fixed_solve(
					data.solver_data,
					B, MatrixX3d(), Beq,
					q);
			}
			else {
				min_quad_with_fixed_solve(
					data.solver_data,
					B, bc, Beq,
					q);
			}


			if (data.energy == RBC_ENERGY_TYPE_RBC) {
				U = data.B * q;
			}
			else {
				U = q;
			}

			if (data.energy == RBC_ENERGY_TYPE_RBC &&
				data.bone_constraint == CONSTRAINED_BONE_CONSTRAINT) {
				fit_hinged_rigid_motion(data.V, data.N, P, I, U);
			}

			// constraint the motion of the bone to be rigid.
			if (data.energy == RBC_ENERGY_TYPE_RBC &&
				(data.bone_constraint == RIGID_BONE_CONSTRAINT ||
				 data.bone_constraint == CONSTRAINED_BONE_CONSTRAINT)) {
				int row = 0;
				for (int i = 1; i < data.N.size(); i++) {
					int nb = data.N(i);
					MatrixXd Ub = U.block(data.nf + row, 0, nb, data.dim);
					MatrixXd Vb = data.V.block(data.nf + row, 0, nb, data.dim);
					Eigen::Matrix3d R;
					Eigen::RowVector3d t;
					Eigen::VectorXd w = Eigen::VectorXd::Ones(nb, 1);
					fit_rigid_motion(Vb, Ub, w, R, t);
					Ub = Vb * R + t.replicate(nb, 1);
					U.block(data.nf + row, 0, nb, data.dim) = Ub;
					row += nb;
				}
			}
			//if (data.collision_enabled) {
			//	double floor_y = data.floor_y;
			//	for (int i = 0; i < U0.rows(); i++) {
			//		if (U0.row(i).y() < floor_y) {
			//			U(i, 1) = floor_y;
			//		}
			//	}
			//}

			iter++;
		}
		if (data.with_dynamics)
		{
			data.vel = (U - U0)/data.h;
		}
		return true;
	}
	
#ifdef IGL_STATIC_LIBRARY
template bool igl::rbc_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, igl::RBCData&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
// generated from msvc error
template bool igl::rbc_solve<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const &, struct igl::RBCData &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);
//
template bool igl::rbc_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, struct igl::RBCData &);
template bool igl::rbc_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, int, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, struct igl::RBCData &);



#endif