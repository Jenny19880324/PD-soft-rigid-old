// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_MIN_QUAD_H
#define IGL_MIN_QUAD_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

namespace igl
{
  template <typename T>
  struct min_quad_data;
  //
  // MIN_QUAD Minimize a quadratic energy of the form
  //
  // trace( 0.5*Z'*A*Z + Z'*B + constant )
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   A  n by n matrix of quadratic coefficients
  //   pd flag specifying whether A(unknown,unknown) is positive definite
  // Outputs:
  //   data  factorization struct with all necessary information to solve
  //     using min_quad_solve
  // Returns true on success, false on error
  //
  //
  template <typename T>
  IGL_INLINE bool min_quad_precompute(
    const Eigen::SparseMatrix<T>& A,
    const bool pd,
    min_quad_with_fixed_data<T> & data
    );
  // Solves a system previously factored using min_quad_with_fixed_precompute
  //
  // Template:
  //   T  type of sparse matrix (e.g. double)
  // Inputs:
  //   data  factorization struct with all necessary precomputation to solve
  //   B  n by k column of linear coefficients
  // Outputs:
  //   Z  n by k solution
  // Returns true on success, false on error
  template <
    typename T,
    typename DerivedB,
    typename DerivedZ>
  IGL_INLINE bool min_quad_solve(
    const min_quad_data<T> & data,
    const Eigen::MatrixBase<DerivedB> & B,
    Eigen::PlainObjectBase<DerivedZ> & Z);

}

template <typename T>
struct igl::min_quad_data
{
  // Size of original system: number of unknowns + number of knowns
  int n;
  // Whether A(unknown,unknown) is positive definite
  bool Auu_pd;
  // Whether A(unknown,unknown) is symmetric
  bool Auu_sym;
  enum SolverType
  {
    LLT = 0,
    LDLT = 1,
    LU = 2,
    QR_LLT = 3,
    NUM_SOLVER_TYPES = 4
  } solver_type;
  // Solvers
  Eigen::SimplicialLLT <Eigen::SparseMatrix<T > > llt;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<T > > ldlt;
  Eigen::SparseLU<Eigen::SparseMatrix<T, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> >   lu;
  // Debug
  Eigen::SparseMatrix<T> NA;
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> NB;
};

#ifndef IGL_STATIC_LIBRARY
#  include "min_quad.cpp"
#endif

#endif
