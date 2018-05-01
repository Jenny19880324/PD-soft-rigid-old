#include "rbc.h"
#include "shape_matrix.h"
#include "massmatrix.h"
#include "fit_rotations.h"
#include "cotmatrix.h"
#include "rbc_rhs.h"

template<
	typename DerivedV,
	typename DerivedF,
	typename Derivedb>
IGL_INLINE bool igl::rbc_precomputation(
	const Eigen::PlainObjectBase<DerivedV> & V,
	const Eigen::PlainObjectBase<DerivedV> & Vb,
	const Eigen::PlainObjectBase<DerivedF> & F,
	const int dim,
	const Eigen::PlainObjectBase<Derivedb> & b,
	RBCData & data)
{
	using namespace std;
	using namespace Eigen;
	typedef typename DerivedV::Scalar Scalar;
	// number of vertices
	const int n = V.rows();
	const int nb = Vb.rows();
	const int nf = n - nb;

	data.n = n;
	data.nf = nf;
	data.nb = nb;
	assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds");
	assert((b.size() == 0 || b.minCoeff() >=0) && "b out of bounds");
	// remember b
	data.b = b;
	data.F = F;
	// dimension
	assert((dim == 3 || dim == 2) && "dim should be 2 or 3");
	data.dim = dim;
	// Defaults
	data.f_ext = MatrixXd::Zero(n, data.dim);
	
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
	SparseMatrix<double> Q = L.eval();

	if (data.with_dynamics)
	{
		const double h = data.h;
		assert(h != 0);
		SparseMatrix<double> M;
		massmatrix(V, F, MASSMATRIX_TYPE_DEFAULT, data.M);
		const double dw = 1. / data.mu;
		SparseMatrix<double> DQ = dw * 1./(h * h) * data.M;
		Q += DQ;
		// Dummy external forces
		data.f_ext = MatrixXd::Zero(n, data.dim);
		data.vel = MatrixXd::Zero(n, data.dim);
	}

	if (eff_energy == RBC_ENERGY_TYPE_RBC) {
		data.Ab.resize(nb, 4 * data.m);
		data.Ab << Vb, MatrixXd::Constant(nb, 1, 1);
		data.B.resize(nf + nb, nf + 4 * data.m);
		std::vector<Triplet<double>> B_IJV;
		for (int i = 0; i < nf; i++) {
			B_IJV.push_back(Triplet<double>(i, i, 1));
		}
		for (int i = 0; i < nb; i++) {
			for (int j = 0; j < 4; j++) {
				B_IJV.push_back(Triplet<double>(nf + i, nf + j, data.Ab(i, j)));
			}
		}
		data.B.setFromTriplets(B_IJV.begin(), B_IJV.end());
		data.B_trans = data.B.transpose();
		Q = data.B_trans * Q * data.B;
		
		data.T.resize(4 * data.m, 3);
		data.T.setZero();
		data.T.block(0, 0, 3, 3) = MatrixXd::Identity(3, 3);
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
				Matrix3d Dm = data.SM.block(3 * r, 0, 3, 3);
				Matrix3d Ds = S.block(3 * r, 0, 3, 3);
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
				Dl = dw * 1./(h*h)*data.M*(-U0-h*data.vel) - data.f_ext;
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

			if (data.energy == RBC_ENERGY_TYPE_RBC) {
				B = data.B_trans * B;
			}

			min_quad_with_fixed_solve(
			data.solver_data,
			B, bc, Beq,
			q);

			if (data.energy == RBC_ENERGY_TYPE_RBC) {
				U = data.B * q;
			}
			else {
				U = q;
			}

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
template bool igl::rbc_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, struct igl::RBCData &);

#endif