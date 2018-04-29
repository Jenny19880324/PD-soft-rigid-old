#ifndef IGL_RBC_H
#define IGL_RBC_H
#include "igl_inline.h"
#include "min_quad_with_fixed.h"
#include "RBCEnergyType.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl
{
	struct RBCData
	{
		// n  #V
		// nf #flesh V
		// nb #bone V
		// m  the number of independent rigid body
		// G  #V list of group indices (1 to k) for each vertex, such that vertex i
		//    is assigned to group G(i)
		// energy type of energy to use
		// with_dynamics whether using dynamics (need to call rbc_precomputation 
		// after changing)
		// f_ext #V by dim list of external forces
		// vel #V by dim list of velocities
		// B #V by (nf + 4m) 
		// T 4m by 3
		// h dynamics time step
		// mu mu = k / (2 * (1 + v)) k is Young's modulus, v is poisson's ratio
		// max_iter maximum inner iterations
		// J rhs pre-multiplier
		// M mass matrix
		// solver_data quadratic solver data
		// b list of boundary indices into V
		// dim dimension being used for solving
		int n, nf, nb, m;
		Eigen::VectorXi G;
		Eigen::MatrixXi F;
		RBCEnergyType energy;
		bool with_dynamics;
		Eigen::MatrixXd f_ext, vel, Ab, T;
		double h;
		double mu;
		int max_iter;
		Eigen::SparseMatrix<double> J, M;
		Eigen::SparseMatrix<double> L;
		Eigen::SparseMatrix<double> B, B_trans;
		Eigen::MatrixXd SM; // restpose shape matrix
		min_quad_with_fixed_data<double> solver_data;
		Eigen::VectorXi b;
		int dim;
		RBCData():
		n(0),
		nf(0),
		nb(0),
		m(1),
		G(),
		energy(RBC_ENERGY_TYPE_DEFAULT),
		with_dynamics(false),
		f_ext(),
		h(1),
		mu(1),
		max_iter(10),
		J(),
		L(),
		solver_data(),
		b(),
		dim(-1) // force this to be set by _precomputation
		{
		};	
	};
	
	// Compute necessary information to start using an ARAP deformation
	//
	// Inputs:
	// V #V by dim list of mesh positions
	// F #F by simplex_size list of triangle|tet indices into V
	// dim dimension being used at solve time. For deformation usually dim =
	// V.cols(), for surface parameterization V.cols() = 3 and dim = 2
	// b #b list of "boundary" fixed vertex indices into V
	// Outputs:
	// data struct containing necessary precomputation
	template <
		typename DerivedV,
		typename DerivedF,
		typename Derivedb>
	IGL_INLINE bool rbc_precomputation(
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedV> & Vb,
		const Eigen::PlainObjectBase<DerivedF> & F,
		const int dim,
		const Eigen::PlainObjectBase<Derivedb> & b,
		RBCData & data);
		
	// Inputs:
	// bc #b by dim list of boundary conditions
	// data struct containing necessary precomputation and parameters
	// U #V by dim initial guess
	template <
		typename Derivedbc,
		typename DerivedU>
	IGL_INLINE bool rbc_solve(
		const Eigen::PlainObjectBase<Derivedbc> &bc,
		RBCData & data,
		Eigen::PlainObjectBase<DerivedU> & U);
}; // namespace igl

#endif