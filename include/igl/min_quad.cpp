// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "min_quad.h"

#include "slice.h"
#include "is_symmetric.h"
#include "find.h"
#include "sparse.h"
#include "repmat.h"
#include "matlab_format.h"
#include "EPS.h"
#include "cat.h"

//#include <Eigen/SparseExtra>
// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>
#include <cassert>
#include <cstdio>
#include <iostream>

template <typename T>
IGL_INLINE bool igl::min_quad_precompute(
  const Eigen::SparseMatrix<T>& A,
  const bool pd,
  min_quad_data<T> & data
  )
{
//#define MIN_QUAD_WITH_FIXED_CPP_DEBUG
  using namespace Eigen;
  using namespace std;
  // number of rows
  int n = A.rows();
  // cache problem size
  data.n = n;

  assert(A.rows() == n && "A should be square");
  assert(A.cols() == n && "A should be square");

  data.llt.compute(A);
  switch(data.llt.info())
  {
	case Eigen::Success:
	  break;
	case Eigen::NumericalIssue:
	  cerr<<"Error: Numerical issue."<<endl;
	  return false;
	default:
	  cerr<<"Error: Other."<<endl;
	  return false;
  }
  data.solver_type = min_quad_with_fixed_data<T>::LLT;
  
  return true;
}


template <
  typename T,
  typename DerivedB,
  typename DerivedY,
  typename DerivedBeq,
  typename DerivedZ,
  typename Derivedsol>
IGL_INLINE bool igl::min_quad_with_fixed_solve(
  const min_quad_with_fixed_data<T> & data,
  const Eigen::MatrixBase<DerivedB> & B,
  const Eigen::MatrixBase<DerivedY> & Y,
  const Eigen::MatrixBase<DerivedBeq> & Beq,
  Eigen::PlainObjectBase<DerivedZ> & Z,
  Eigen::PlainObjectBase<Derivedsol> & sol)
{
  using namespace std;
  using namespace Eigen;
  typedef Matrix<T,Dynamic,1> VectorXT;
  typedef Matrix<T,Dynamic,Dynamic> MatrixXT;
  // number of known rows
  int kr = data.known.size();
  if(kr!=0)
  {
    assert(kr == Y.rows());
  }
  // number of columns to solve
  int cols = Y.cols();
  assert(B.cols() == 1 || B.cols() == cols);
  assert(Beq.size() == 0 || Beq.cols() == 1 || Beq.cols() == cols);

  // resize output
  Z.resize(data.n,cols);
  // Set known values
  for(int i = 0;i < kr;i++)
  {
    for(int j = 0;j < cols;j++)
    {
      Z(data.known(i),j) = Y(i,j);
    }
  }

  if(data.Aeq_li)
  {
    // number of lagrange multipliers aka linear equality constraints
    int neq = data.lagrange.size();
    // append lagrange multiplier rhs's
    MatrixXT BBeq(B.rows() + Beq.rows(),cols);
    if(B.size() > 0)
    {
      BBeq.topLeftCorner(B.rows(),cols) = B.replicate(1,B.cols()==cols?1:cols);
    }
    if(Beq.size() > 0)
    {
      BBeq.bottomLeftCorner(Beq.rows(),cols) = -2.0*Beq.replicate(1,Beq.cols()==cols?1:cols);
    }

    // Build right hand side
    MatrixXT BBequlcols;
    igl::slice(BBeq,data.unknown_lagrange,1,BBequlcols);
    MatrixXT NB;
    if(kr == 0)
    {
      NB = BBequlcols;
    }else
    {
      NB = data.preY * Y + BBequlcols;
    }

    //std::cout<<"NB=["<<std::endl<<NB<<std::endl<<"];"<<std::endl;
    //cout<<matlab_format(NB,"NB")<<endl;
    switch(data.solver_type)
    {
      case igl::min_quad_with_fixed_data<T>::LLT:
        sol = data.llt.solve(NB);
        break;
      case igl::min_quad_with_fixed_data<T>::LDLT:
        sol = data.ldlt.solve(NB);
        break;
      case igl::min_quad_with_fixed_data<T>::LU:
        // Not a bottleneck
        sol = data.lu.solve(NB);
        break;
      default:
        cerr<<"Error: invalid solver type"<<endl;
        return false;
    }
    //std::cout<<"sol=["<<std::endl<<sol<<std::endl<<"];"<<std::endl;
    // Now sol contains sol/-0.5
    sol *= -0.5;
    // Now sol contains solution
    // Place solution in Z
    for(int i = 0;i<(sol.rows()-neq);i++)
    {
      for(int j = 0;j<sol.cols();j++)
      {
        Z(data.unknown_lagrange(i),j) = sol(i,j);
      }
    }
  }else
  {
    assert(data.solver_type == min_quad_with_fixed_data<T>::QR_LLT);
    MatrixXT eff_Beq;
    // Adjust Aeq rhs to include known parts
    eff_Beq =
      //data.AeqTQR.colsPermutation().transpose() * (-data.Aeqk * Y + Beq);
      data.AeqTET * (-data.Aeqk * Y + Beq.replicate(1,Beq.cols()==cols?1:cols));
    // Where did this -0.5 come from? Probably the same place as above.
    MatrixXT Bu;
    slice(B,data.unknown,1,Bu);
    MatrixXT NB;
    NB = -0.5*(Bu.replicate(1,B.cols()==cols?1:cols) + data.preY * Y);
    // Trim eff_Beq
    const int nc = data.AeqTQR.rank();
    const int neq = Beq.rows();
    eff_Beq = eff_Beq.topLeftCorner(nc,cols).eval();
    data.AeqTR1T.template triangularView<Lower>().solveInPlace(eff_Beq);
    // Now eff_Beq = (data.AeqTR1T \ (data.AeqTET * (-data.Aeqk * Y + Beq)))
    MatrixXT lambda_0;
    lambda_0 = data.AeqTQ1 * eff_Beq;
    //cout<<matlab_format(lambda_0,"lambda_0")<<endl;
    MatrixXT QRB;
    QRB = -data.AeqTQ2T * (data.Auu * lambda_0) + data.AeqTQ2T * NB;
    Derivedsol lambda;
    lambda = data.llt.solve(QRB);
    // prepare output
    Derivedsol solu;
    solu = data.AeqTQ2 * lambda + lambda_0;
    //  http://www.math.uh.edu/~rohop/fall_06/Chapter3.pdf
    Derivedsol solLambda;
    {
      Derivedsol temp1,temp2;
      temp1 = (data.AeqTQ1T * NB - data.AeqTQ1T * data.Auu * solu);
      data.AeqTR1.template triangularView<Upper>().solveInPlace(temp1);
      //cout<<matlab_format(temp1,"temp1")<<endl;
      temp2 = Derivedsol::Zero(neq,cols);
      temp2.topLeftCorner(nc,cols) = temp1;
      //solLambda = data.AeqTQR.colsPermutation() * temp2;
      solLambda = data.AeqTE * temp2;
    }
    // sol is [Z(unknown);Lambda]
    assert(data.unknown.size() == solu.rows());
    assert(cols == solu.cols());
    assert(data.neq == neq);
    assert(data.neq == solLambda.rows());
    assert(cols == solLambda.cols());
    sol.resize(data.unknown.size()+data.neq,cols);
    sol.block(0,0,solu.rows(),solu.cols()) = solu;
    sol.block(solu.rows(),0,solLambda.rows(),solLambda.cols()) = solLambda;
    for(int u = 0;u<data.unknown.size();u++)
    {
      for(int j = 0;j<Z.cols();j++)
      {
        Z(data.unknown(u),j) = solu(u,j);
      }
    }
  }
  return true;
}

template <
  typename T,
  typename DerivedB,
  typename DerivedY,
  typename DerivedBeq,
  typename DerivedZ>
IGL_INLINE bool igl::min_quad_with_fixed_solve(
  const min_quad_with_fixed_data<T> & data,
  const Eigen::MatrixBase<DerivedB> & B,
  const Eigen::MatrixBase<DerivedY> & Y,
  const Eigen::MatrixBase<DerivedBeq> & Beq,
  Eigen::PlainObjectBase<DerivedZ> & Z)
{
  Eigen::Matrix<typename DerivedZ::Scalar, Eigen::Dynamic, Eigen::Dynamic> sol;
  return min_quad_with_fixed_solve(data,B,Y,Beq,Z,sol);
}

template <
  typename T,
  typename Derivedknown,
  typename DerivedB,
  typename DerivedY,
  typename DerivedBeq,
  typename DerivedZ>
IGL_INLINE bool igl::min_quad_with_fixed(
  const Eigen::SparseMatrix<T>& A,
  const Eigen::MatrixBase<DerivedB> & B,
  const Eigen::MatrixBase<Derivedknown> & known,
  const Eigen::MatrixBase<DerivedY> & Y,
  const Eigen::SparseMatrix<T>& Aeq,
  const Eigen::MatrixBase<DerivedBeq> & Beq,
  const bool pd,
  Eigen::PlainObjectBase<DerivedZ> & Z)
{
  min_quad_with_fixed_data<T> data;
  if(!min_quad_with_fixed_precompute(A,known,Aeq,pd,data))
  {
    return false;
  }
  return min_quad_with_fixed_solve(data,B,Y,Beq,Z);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template bool igl::min_quad_with_fixed<double, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template bool igl::min_quad_with_fixed_precompute<double, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int> const&, bool, igl::min_quad_with_fixed_data<double>&);
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::min_quad_with_fixed_precompute<double, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int> const&, bool, igl::min_quad_with_fixed_data<double>&);
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::min_quad_with_fixed<double, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
// generated from msvc error
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(struct igl::min_quad_with_fixed_data<double> const &, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const &, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(struct igl::min_quad_with_fixed_data<double> const &, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const &, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const &, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);
#endif
