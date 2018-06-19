// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_READMESH_H
#define IGL_READMESH_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <unordered_map>

namespace igl
{
  // load a tetrahedral volume mesh from a .mesh file
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Input:
  //   mesh_file_name  path of .mesh file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   T  #T list of tet indices into vertex positions
  //   F  #F list of face indices into vertex positions
  //
  // Known bugs: Holes and regions are not supported
  template <typename Scalar, typename Index>
  IGL_INLINE bool readMESH(
    const std::string mesh_file_name,
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Index > > & T,
    std::vector<std::vector<Index > > & F);
  // Inputs:
  //   mesh_file  pointer to already opened .mesh file 
  // Outputs:
  //   mesh_file  closed file
  template <typename Scalar, typename Index>
  IGL_INLINE bool readMESH(
    FILE * mesh_file,
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Index > > & T,
    std::vector<std::vector<Index > > & F);

  // Input:
  //   mesh_file_name  path of .mesh file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   T  eigen int matrix #T by 4
  //   F  eigen int matrix #F by 3
  template <typename DerivedV, typename DerivedF, typename DerivedT>
  IGL_INLINE bool readMESH(
    const std::string mesh_file_name,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedF>& F);
  // Inputs:
  //   mesh_file  pointer to already opened .mesh file 
  // Outputs:
  //   mesh_file  closed file
  template <typename DerivedV, typename DerivedF, typename DerivedT>
  IGL_INLINE bool readMESH(
    FILE * mesh_file,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedF>& F);


	// Input:
	// mesh_file_name path of .mesh file
	// Outputs:
	// V eigen double matrix #V by 3
	// T eigen int matrix #T by 4
	// F eigen int matrix #F by 3
	// C eigen double matrix #F by 4
    // N eigen int matrix #region by 1 
    // (number_of_elastic_vertices, number_of_first_rigid_body_vertices, number_of_second_rigid_body_vertices ...)
	template <
		typename DerivedV, 
		typename DerivedT, 
		typename DerivedF, 
		typename DerivedC,
		typename DerivedN>
	IGL_INLINE bool readMESH(
		const std::string mesh_file_name,
		Eigen::PlainObjectBase<DerivedV>& V,
		Eigen::PlainObjectBase<DerivedT>& T,
		Eigen::PlainObjectBase<DerivedF>& F,
		Eigen::PlainObjectBase<DerivedC>& C,
		Eigen::PlainObjectBase<DerivedN>& N,
		Eigen::VectorXi & A);
	// Inputs:
	// mesh_file pointer to already opened .mesh file
	// Outputs:
	// mesh_file closed file
	template <
		typename DerivedV, 
		typename DerivedT, 
		typename DerivedF, 
		typename DerivedC,
		typename DerivedN>
	IGL_INLINE bool readMESH(
		FILE * mesh_file,
		Eigen::PlainObjectBase<DerivedV> & V,
		Eigen::PlainObjectBase<DerivedT> & T,
		Eigen::PlainObjectBase<DerivedF> & F,
		Eigen::PlainObjectBase<DerivedC> & C,
		Eigen::PlainObjectBase<DerivedN> & N,
		Eigen::VectorXi & A);

	// readMESH for self collision
	// Outputs:
	// FF, triangle faces without duplicate.
	template <
		typename DerivedV,
		typename DerivedT,
		typename DerivedF,
		typename DerivedC,
		typename DerivedN>
		IGL_INLINE bool readMESH(
			const std::string mesh_file_name,
			Eigen::PlainObjectBase< DerivedV > & V,
			Eigen::PlainObjectBase< DerivedT > & T,
			Eigen::PlainObjectBase< DerivedF > & F,
			Eigen::PlainObjectBase< DerivedF > & SF,
			Eigen::PlainObjectBase< DerivedC > & C,
			Eigen::PlainObjectBase< DerivedN > & N,
			Eigen::VectorXi & A);

	template <
		typename DerivedV,
		typename DerivedT,
		typename DerivedF,
		typename DerivedC,
		typename DerivedN>
		IGL_INLINE bool readMESH(
			FILE * mesh_file,
			Eigen::PlainObjectBase< DerivedV > & V,
			Eigen::PlainObjectBase< DerivedT > & T,
			Eigen::PlainObjectBase< DerivedF > & F,
			Eigen::PlainObjectBase< DerivedF > & SF,
			Eigen::PlainObjectBase< DerivedC > & C,
			Eigen::PlainObjectBase< DerivedN > & N,
			Eigen::VectorXi & A);


	// readMESH for self collision
	// Outputs:
	// SF surface facets 
	// SV surface vertices
	template <
		typename DerivedV,
		typename DerivedT,
		typename DerivedF,
		typename DerivedC,
		typename DerivedN>
		IGL_INLINE bool readMESH(
			const std::string mesh_file_name,
			Eigen::PlainObjectBase<DerivedV> & V,
			Eigen::PlainObjectBase<DerivedT> & T,
			Eigen::PlainObjectBase<DerivedF> & F,
			Eigen::PlainObjectBase<DerivedT> & ST,
			Eigen::PlainObjectBase<DerivedF> & SF,
			Eigen::PlainObjectBase<DerivedC> & C,
			Eigen::PlainObjectBase<DerivedN> & N,
			Eigen::VectorXi & A,
			Eigen::VectorXi & SV);

	template <
		typename DerivedV,
		typename DerivedT,
		typename DerivedF,
		typename DerivedC,
		typename DerivedN>
		IGL_INLINE bool readMESH(
			FILE * mesh_file,
			Eigen::PlainObjectBase<DerivedV> & V,
			Eigen::PlainObjectBase<DerivedT> & T,
			Eigen::PlainObjectBase<DerivedF> & F,
			Eigen::PlainObjectBase<DerivedT> & ST,
			Eigen::PlainObjectBase<DerivedF> & SF,
			Eigen::PlainObjectBase<DerivedC> & C,
			Eigen::PlainObjectBase<DerivedN> & N,
			Eigen::VectorXi & A,
			Eigen::VectorXi & SV);

}

#ifndef IGL_STATIC_LIBRARY
#  include "readMESH.cpp"
#endif

#endif
