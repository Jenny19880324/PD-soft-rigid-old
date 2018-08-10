// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "readMESH.h"
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <iostream>

template <typename Scalar, typename Index>
IGL_INLINE bool igl::readMESH(
  const std::string mesh_file_name,
  std::vector<std::vector<Scalar > > & V,
  std::vector<std::vector<Index > > & T,
  std::vector<std::vector<Index > > & F)
{
  using namespace std;
  FILE * mesh_file = fopen(mesh_file_name.c_str(),"r");
  if(NULL==mesh_file)
  {
    fprintf(stderr,"IOError: %s could not be opened...",mesh_file_name.c_str());
    return false;
  }
  return igl::readMESH(mesh_file,V,T,F);
}

template <typename Scalar, typename Index>
IGL_INLINE bool igl::readMESH(
  FILE * mesh_file,
  std::vector<std::vector<Scalar > > & V,
  std::vector<std::vector<Index > > & T,
  std::vector<std::vector<Index > > & F)
{
  using namespace std;
#ifndef LINE_MAX
#  define LINE_MAX 2048
#endif
  char line[LINE_MAX];
  bool still_comments;
  V.clear();
  T.clear();
  F.clear();

  // eat comments at beginning of file
  still_comments= true;
  while(still_comments)
  {
    if(fgets(line,LINE_MAX,mesh_file) == NULL)
    {
      fprintf(stderr, "Error: couldn't find start of .mesh file");
      fclose(mesh_file);
      return false;
    }
    still_comments = (line[0] == '#' || line[0] == '\n');
  }
  char str[LINE_MAX];
  sscanf(line," %s",str);
  // check that first word is MeshVersionFormatted
  if(0!=strcmp(str,"MeshVersionFormatted"))
  {
    fprintf(stderr,
      "Error: first word should be MeshVersionFormatted not %s\n",str);
    fclose(mesh_file);
    return false;
  }

  int one = -1;
  if(2 != sscanf(line,"%s %d",str,&one))
  {
    // 1 appears on next line?
    fscanf(mesh_file," %d",&one);
  }
  if(one != 1)
  {
    fprintf(stderr,"Error: second word should be 1 not %d\n",one);
    fclose(mesh_file);
    return false;
  }

  // eat comments
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  sscanf(line," %s",str);
  // check that third word is Dimension
  if(0!=strcmp(str,"Dimension"))
  {
    fprintf(stderr,"Error: third word should be Dimension not %s\n",str);
    fclose(mesh_file);
    return false;
  }
  int three = -1;
  if(2 != sscanf(line,"%s %d",str,&three))
  {
    // 1 appears on next line?
    fscanf(mesh_file," %d",&three);
  }
  if(three != 3)
  {
    fprintf(stderr,"Error: only Dimension 3 supported not %d\n",three);
    fclose(mesh_file);
    return false;
  }

  // eat comments
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  sscanf(line," %s",str);
  // check that fifth word is Vertices
  if(0!=strcmp(str,"Vertices"))
  {
    fprintf(stderr,"Error: fifth word should be Vertices not %s\n",str);
    fclose(mesh_file);
    return false;
  }

  //fgets(line,LINE_MAX,mesh_file);

  int number_of_vertices;
  if(1 != fscanf(mesh_file," %d",&number_of_vertices) || number_of_vertices > 1000000000)
  {
    fprintf(stderr,"Error: expecting number of vertices less than 10^9...\n");
    fclose(mesh_file);
    return false;
  }
  // allocate space for vertices
  V.resize(number_of_vertices,vector<Scalar>(3,0));
  int extra;
  for(int i = 0;i<number_of_vertices;i++)
  {
    double x,y,z;
    if(4 != fscanf(mesh_file," %lg %lg %lg %d",&x,&y,&z,&extra))
    {
      fprintf(stderr,"Error: expecting vertex position...\n");
      fclose(mesh_file);
      return false;
    }
    V[i][0] = x;
    V[i][1] = y;
    V[i][2] = z;
  }

  // eat comments
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  sscanf(line," %s",str);
  // check that sixth word is Triangles
  if(0!=strcmp(str,"Triangles"))
  {
    fprintf(stderr,"Error: sixth word should be Triangles not %s\n",str);
    fclose(mesh_file);
    return false;
  }
  int number_of_triangles;
  if(1 != fscanf(mesh_file," %d",&number_of_triangles))
  {
    fprintf(stderr,"Error: expecting number of triangles...\n");
    fclose(mesh_file);
    return false;
  }
  // allocate space for triangles
  F.resize(number_of_triangles,vector<Index>(3));
  // triangle indices
  int tri[3];
  for(int i = 0;i<number_of_triangles;i++)
  {
    if(4 != fscanf(mesh_file," %d %d %d %d",&tri[0],&tri[1],&tri[2],&extra))
    {
      printf("Error: expecting triangle indices...\n");
      return false;
    }
    for(int j = 0;j<3;j++)
    {
      F[i][j] = tri[j]-1;
    }
  }

  // eat comments
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  sscanf(line," %s",str);
  // check that sixth word is Triangles
  if(0!=strcmp(str,"Tetrahedra"))
  {
    fprintf(stderr,"Error: seventh word should be Tetrahedra not %s\n",str);
    fclose(mesh_file);
    return false;
  }
  int number_of_tetrahedra;
  if(1 != fscanf(mesh_file," %d",&number_of_tetrahedra))
  {
    fprintf(stderr,"Error: expecting number of tetrahedra...\n");
    fclose(mesh_file);
    return false;
  }
  // allocate space for tetrahedra
  T.resize(number_of_tetrahedra,vector<Index>(4));
  // tet indices
  int a,b,c,d;
  for(int i = 0;i<number_of_tetrahedra;i++)
  {
    if(5 != fscanf(mesh_file," %d %d %d %d %d",&a,&b,&c,&d,&extra))
    {
      fprintf(stderr,"Error: expecting tetrahedra indices...\n");
      fclose(mesh_file);
      return false;
    }
    T[i][0] = a-1;
    T[i][1] = b-1;
    T[i][2] = c-1;
    T[i][3] = d-1;
  }
  fclose(mesh_file);
  return true;
}

#include <Eigen/Core>
#include "list_to_matrix.h"


template <typename DerivedV, typename DerivedF, typename DerivedT>
IGL_INLINE bool igl::readMESH(
  const std::string mesh_file_name,
  Eigen::PlainObjectBase<DerivedV>& V,
  Eigen::PlainObjectBase<DerivedT>& T,
  Eigen::PlainObjectBase<DerivedF>& F)
{
  using namespace std;
  FILE * mesh_file = fopen(mesh_file_name.c_str(),"r");
  if(NULL==mesh_file)
  {
    fprintf(stderr,"IOError: %s could not be opened...",mesh_file_name.c_str());
    return false;
  }
  return readMESH(mesh_file,V,T,F);
}

template <typename DerivedV, typename DerivedF, typename DerivedT>
IGL_INLINE bool igl::readMESH(
  FILE * mesh_file,
  Eigen::PlainObjectBase<DerivedV>& V,
  Eigen::PlainObjectBase<DerivedT>& T,
  Eigen::PlainObjectBase<DerivedF>& F)
{
  using namespace std;
#ifndef LINE_MAX
#  define LINE_MAX 2048
#endif
  char line[LINE_MAX];
  bool still_comments;

  // eat comments at beginning of file
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  char str[LINE_MAX];
  sscanf(line," %s",str);
  // check that first word is MeshVersionFormatted
  if(0!=strcmp(str,"MeshVersionFormatted"))
  {
    fprintf(stderr,
      "Error: first word should be MeshVersionFormatted not %s\n",str);
    fclose(mesh_file);
    return false;
  }
  int one = -1;
  if(2 != sscanf(line,"%s %d",str,&one))
  {
    // 1 appears on next line?
    fscanf(mesh_file," %d",&one);
  }
  if(one != 1)
  {
    fprintf(stderr,"Error: second word should be 1 not %d\n",one);
    fclose(mesh_file);
    return false;
  }

  // eat comments
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  sscanf(line," %s",str);
  // check that third word is Dimension
  if(0!=strcmp(str,"Dimension"))
  {
    fprintf(stderr,"Error: third word should be Dimension not %s\n",str);
    fclose(mesh_file);
    return false;
  }
  int three = -1;
  if(2 != sscanf(line,"%s %d",str,&three))
  {
    // 1 appears on next line?
    fscanf(mesh_file," %d",&three);
  }
  if(three != 3)
  {
    fprintf(stderr,"Error: only Dimension 3 supported not %d\n",three);
    fclose(mesh_file);
    return false;
  }

  // eat comments
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  sscanf(line," %s",str);
  // check that fifth word is Vertices
  if(0!=strcmp(str,"Vertices"))
  {
    fprintf(stderr,"Error: fifth word should be Vertices not %s\n",str);
    fclose(mesh_file);
    return false;
  }

  //fgets(line,LINE_MAX,mesh_file);

  int number_of_vertices;
  if(1 != fscanf(mesh_file," %d",&number_of_vertices) || number_of_vertices > 1000000000)
  {
    fprintf(stderr,"Error: expecting number of vertices less than 10^9...\n");
    fclose(mesh_file);
    return false;
  }
  // allocate space for vertices
  V.resize(number_of_vertices,3);
  int extra;
  for(int i = 0;i<number_of_vertices;i++)
  {
    double x,y,z;
    if(4 != fscanf(mesh_file," %lg %lg %lg %d",&x,&y,&z,&extra))
    {
      fprintf(stderr,"Error: expecting vertex position...\n");
      fclose(mesh_file);
      return false;
    }
    V(i,0) = x;
    V(i,1) = y;
    V(i,2) = z;
  }

  // eat comments
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  sscanf(line," %s",str);
  // check that sixth word is Triangles
  if(0!=strcmp(str,"Triangles"))
  {
    fprintf(stderr,"Error: sixth word should be Triangles not %s\n",str);
    fclose(mesh_file);
    return false;
  }
  int number_of_triangles;
  if(1 != fscanf(mesh_file," %d",&number_of_triangles))
  {
    fprintf(stderr,"Error: expecting number of triangles...\n");
    fclose(mesh_file);
    return false;
  }
  // allocate space for triangles
  F.resize(number_of_triangles,3);
  // triangle indices
  int tri[3];
  for(int i = 0;i<number_of_triangles;i++)
  {
    if(4 != fscanf(mesh_file," %d %d %d %d",&tri[0],&tri[1],&tri[2],&extra))
    {
      printf("Error: expecting triangle indices...\n");
      return false;
    }
    for(int j = 0;j<3;j++)
    {
      F(i,j) = tri[j]-1;
    }
  }

  // eat comments
  still_comments= true;
  while(still_comments)
  {
    fgets(line,LINE_MAX,mesh_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }

  sscanf(line," %s",str);
  // check that sixth word is Triangles
  if(0!=strcmp(str,"Tetrahedra"))
  {
    fprintf(stderr,"Error: seventh word should be Tetrahedra not %s\n",str);
    fclose(mesh_file);
    return false;
  }
  int number_of_tetrahedra;
  if(1 != fscanf(mesh_file," %d",&number_of_tetrahedra))
  {
    fprintf(stderr,"Error: expecting number of tetrahedra...\n");
    fclose(mesh_file);
    return false;
  }
  // allocate space for tetrahedra
  T.resize(number_of_tetrahedra,4);
  // tet indices
  int a,b,c,d;
  for(int i = 0;i<number_of_tetrahedra;i++)
  {
    if(5 != fscanf(mesh_file," %d %d %d %d %d",&a,&b,&c,&d,&extra))
    {
      fprintf(stderr,"Error: expecting tetrahedra indices...\n");
      fclose(mesh_file);
      return false;
    }
    T(i,0) = a-1;
    T(i,1) = b-1;
    T(i,2) = c-1;
    T(i,3) = d-1;
  }
  fclose(mesh_file);
  return true;
}
//{
//  std::vector<std::vector<double> > vV,vT,vF;
//  bool success = igl::readMESH(mesh_file_name,vV,vT,vF);
//  if(!success)
//  {
//    // readMESH already printed error message to std err
//    return false;
//  }
//  bool V_rect = igl::list_to_matrix(vV,V);
//  if(!V_rect)
//  {
//    // igl::list_to_matrix(vV,V) already printed error message to std err
//    return false;
//  }
//  bool T_rect = igl::list_to_matrix(vT,T);
//  if(!T_rect)
//  {
//    // igl::list_to_matrix(vT,T) already printed error message to std err
//    return false;
//  }
//  bool F_rect = igl::list_to_matrix(vF,F);
//  if(!F_rect)
//  {
//    // igl::list_to_matrix(vF,F) already printed error message to std err
//    return false;
//  }
//  assert(V.cols() == 3);
//  assert(T.cols() == 4);
//  assert(F.cols() == 3);
//  return true;
//}

template <
	typename DerivedV, 
	typename DerivedT, 
	typename DerivedF, 
	typename DerivedC,
	typename DerivedN>
IGL_INLINE bool igl::readMESH(
	const std::string mesh_file_name,
	Eigen::PlainObjectBase<DerivedV>& V,
	Eigen::PlainObjectBase<DerivedT>& T,
	Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedC>& C,
	Eigen::PlainObjectBase<DerivedN>& N,
	Eigen::VectorXi & A)
{
	using namespace std;
	FILE * mesh_file = fopen(mesh_file_name.c_str(), "r");
	if (NULL == mesh_file)
	{
		fprintf(stderr, "IOError: %s could not be opened...", mesh_file_name.c_str());
		return false;
	}
	return readMESH(mesh_file, V, T, F, C, N, A);
}


template <
	typename DerivedV, 
	typename DerivedT, 
	typename DerivedF, 
	typename DerivedC,
	typename DerivedN>
IGL_INLINE bool igl::readMESH(
	FILE * mesh_file,
	Eigen::PlainObjectBase<DerivedV>& V,
	Eigen::PlainObjectBase<DerivedT>& T,
	Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedC>& C,
	Eigen::PlainObjectBase<DerivedN>& N,
	Eigen::VectorXi & A)
{
	V.resize(0, Eigen::NoChange);
	T.resize(0, Eigen::NoChange);
	F.resize(0, Eigen::NoChange);
	C.resize(0, Eigen::NoChange);
	N.resize(0, Eigen::NoChange);
	A.resize(0, Eigen::NoChange);

	using namespace std;
#ifndef LINE_MAX
#  define LINE_MAX 2048
#endif
	char line[LINE_MAX];
	bool still_comments;

	// eat comments at beginning of file
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	char str[LINE_MAX];
	sscanf(line, " %s", str);
	// check that first word is MeshVersionFormatted
	if (0 != strcmp(str, "MeshVersionFormatted"))
	{
		fprintf(stderr,
			"Error: first word should be MeshVersionFormatted not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int one = -1;
	if (2 != sscanf(line, "%s %d", str, &one))
	{
		// 1 appears on next line?
		fscanf(mesh_file, " %d", &one);
	}
	if (one != 1)
	{
		fprintf(stderr, "Error: second word should be 1 not %d\n", one);
		fclose(mesh_file);
		return false;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that third word is Dimension
	if (0 != strcmp(str, "Dimension"))
	{
		fprintf(stderr, "Error: third word should be Dimension not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int three = -1;
	if (2 != sscanf(line, "%s %d", str, &three))
	{
		// 1 appears on next line?
		fscanf(mesh_file, " %d", &three);
	}
	if (three != 3)
	{
		fprintf(stderr, "Error: only Dimension 3 supported not %d\n", three);
		fclose(mesh_file);
		return false;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that fifth word is Vertices
	if (0 != strcmp(str, "Vertices"))
	{
		fprintf(stderr, "Error: fifth word should be Vertices not %s\n", str);
		fclose(mesh_file);
		return false;
	}

	//fgets(line,LINE_MAX,mesh_file);

	int number_of_vertices;
	if (1 != fscanf(mesh_file, " %d", &number_of_vertices) || number_of_vertices > 1000000000)
	{
		fprintf(stderr, "Error: expecting number of vertices less than 10^9...\n");
		fclose(mesh_file);
		return false;
	}
	// allocate space for vertices
	V.resize(number_of_vertices, 3);
	int extra;
	for (int i = 0; i<number_of_vertices; i++)
	{
		double x, y, z;
		if (4 != fscanf(mesh_file, " %lg %lg %lg %d", &x, &y, &z, &extra))
		{
			fprintf(stderr, "Error: expecting vertex position...\n");
			fclose(mesh_file);
			return false;
		}
		V(i, 0) = x;
		V(i, 1) = y;
		V(i, 2) = z;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that sixth word is Triangles
	if (0 != strcmp(str, "Triangles"))
	{
		fprintf(stderr, "Error: sixth word should be Triangles not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int number_of_triangles;
	if (1 != fscanf(mesh_file, " %d", &number_of_triangles))
	{
		fprintf(stderr, "Error: expecting number of triangles...\n");
		fclose(mesh_file);
		return false;
	}

	// triangle indices
	int tri[3];
	for (int i = 0; i<number_of_triangles; i++)
	{
		if (4 != fscanf(mesh_file, " %d %d %d %d", &tri[0], &tri[1], &tri[2], &extra))
		{
			printf("Error: expecting triangle indices...\n");
			return false;
		}
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that sixth word is Triangles
	if (0 != strcmp(str, "Tetrahedra"))
	{
		fprintf(stderr, "Error: seventh word should be Tetrahedra not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int number_of_tetrahedra;
	if (1 != fscanf(mesh_file, " %d", &number_of_tetrahedra))
	{
		fprintf(stderr, "Error: expecting number of tetrahedra...\n");
		fclose(mesh_file);
		return false;
	}
	// allocate space for tetrahedra
	A.resize(number_of_tetrahedra); // region attribute
	T.resize(number_of_tetrahedra, 4);
	// tet indices
	int a, b, c, d;
	//[region]->[vertex index]
	// region == 0 => elastic tet, flesh tet
	// region < 0  => rigid tet, bone tet
	std::map<int, std::set<int>> rigid_vertex_map; 
	std::unordered_set<int> rigid_vertex_set;
	for (int i = 0; i<number_of_tetrahedra; i++)
	{
		if (5 != fscanf(mesh_file, " %d %d %d %d %d", &a, &b, &c, &d, &extra))
		{
			fprintf(stderr, "Error: expecting tetrahedra indices...\n");
			fclose(mesh_file);
			return false;
		}
		T(i, 0) = a - 1;
		T(i, 1) = b - 1;
		T(i, 2) = c - 1;
		T(i, 3) = d - 1;
		A(i) = extra;

		if (extra < 0 ) {
			rigid_vertex_map[extra].insert(T(i, 0));
			rigid_vertex_map[extra].insert(T(i, 1));
			rigid_vertex_map[extra].insert(T(i, 2));
			rigid_vertex_map[extra].insert(T(i, 3));

			rigid_vertex_set.insert(T(i, 0));
			rigid_vertex_set.insert(T(i, 1));
			rigid_vertex_set.insert(T(i, 2));
			rigid_vertex_set.insert(T(i, 3));
		}
	}
	fclose(mesh_file);

	// rearrage the order of the vertices so that V = [Vf; Vt]
	std::unordered_map<int, int> RI;

	N.resize(rigid_vertex_map.size() + 1, Eigen::NoChange);
	int idx = 1;
	int nf = number_of_vertices;
	for (auto it = rigid_vertex_map.begin(); it != rigid_vertex_map.end(); it++)
	{
		N(idx) = it->second.size();
		nf -= N(idx);
		idx++;
	}
	N(0) = nf;

	Eigen::PlainObjectBase<DerivedV> rearranged_V;
	rearranged_V.resize(number_of_vertices, 3);
	int number_of_vertices_prev_set = 0;
	for (auto m_it = rigid_vertex_map.begin(); m_it != rigid_vertex_map.end(); m_it++)
	{
		if (m_it != rigid_vertex_map.begin()) {
			number_of_vertices_prev_set += std::prev(m_it)->second.size();
		}
		std::set<int> &one_rigid_set = m_it->second;
		for (auto s_it = one_rigid_set.begin(); s_it != one_rigid_set.end(); s_it++)
		{
			int cnt = number_of_vertices_prev_set + std::distance(one_rigid_set.begin(), s_it);
			rearranged_V.row(nf + cnt) = V.row(*s_it);
			RI[*s_it] = nf + cnt;
		}
	}

	int elastic_idx = 0;
	for (int i = 0; i < number_of_vertices; i++) {
		if (rigid_vertex_set.find(i) == rigid_vertex_set.end()) //elastic vertex
		{
			rearranged_V.row(elastic_idx) = V.row(i);
			RI[i] = elastic_idx;
			elastic_idx++;
		}
	}
	assert(elastic_idx == nf);
	V = rearranged_V;

	// update index of tetrahedron because V is rearranged
	for (int i = 0; i < number_of_tetrahedra; i++) {
		for (int j = 0; j < 4; j++) {
			T(i, j) = RI[T(i, j)];
		}
	}

	// construct Triangle from Tetrahedron
	F.resize(4 * number_of_tetrahedra, 3);
	for (int i = 0; i < number_of_tetrahedra; ++i) {
		F(i * 4 + 0, 0) = T(i, 0); F(i * 4 + 0, 1) = T(i, 1); F(i * 4 + 0, 2) = T(i, 3);
		F(i * 4 + 1, 0) = T(i, 1); F(i * 4 + 1, 1) = T(i, 2); F(i * 4 + 1, 2) = T(i, 3);
		F(i * 4 + 2, 0) = T(i, 2); F(i * 4 + 2, 1) = T(i, 0); F(i * 4 + 2, 2) = T(i, 3);
		F(i * 4 + 3, 0) = T(i, 0); F(i * 4 + 3, 1) = T(i, 2); F(i * 4 + 3, 2) = T(i, 1);
	}

	// assign different colors to faces belong to rigid tets

	// http://colorbrewer2.org
	// 12 data classes
	// qualitative color scheme
	Eigen::MatrixX4d colorScheme;
	colorScheme.resize(12, 4);

	colorScheme <<
		166, 206, 227, 255
		, 31, 120, 180, 255
		, 178, 223, 138, 255
		, 51, 160, 44, 255
		, 251, 154, 153, 255
		, 227, 26, 28, 255
		, 253, 191, 111, 255
		, 255, 127, 0, 255
		, 202, 178, 214, 255
		, 106, 61, 154, 255
		, 255, 255, 153, 255
		, 177, 89, 40, 255;

	colorScheme /= 255.0f;
	C.resize(F.rows(), 4);
	for (int f_i = 0; f_i < F.rows(); ++f_i) {
		const int colorIdx = abs(A(f_i / 4)) % 12;
		Eigen::Vector4d color = colorScheme.row(colorIdx);
		C.row(f_i) = color;
	}
	return true;
}


// readMESH for self collision
template <
	typename DerivedV,
	typename DerivedT,
	typename DerivedF,
	typename DerivedC,
	typename DerivedN>
	IGL_INLINE bool igl::readMESH(
		const std::string mesh_file_name,
		Eigen::PlainObjectBase< DerivedV >& V,
		Eigen::PlainObjectBase< DerivedT >& T,
		Eigen::PlainObjectBase< DerivedF >& F,
		Eigen::PlainObjectBase< DerivedF > & SF,
		Eigen::PlainObjectBase< DerivedC >& C,
		Eigen::PlainObjectBase< DerivedN >& N,
		Eigen::VectorXi & A)
{
	using namespace std;
	FILE * mesh_file = fopen(mesh_file_name.c_str(), "r");
	if (NULL == mesh_file)
	{
		fprintf(stderr, "IOError: %s could not be opened...", mesh_file_name.c_str());
		return false;
	}
	return readMESH(mesh_file, V, T, F, SF, C, N, A);
}


template <
	typename DerivedV,
	typename DerivedT,
	typename DerivedF,
	typename DerivedC,
	typename DerivedN>
	IGL_INLINE bool igl::readMESH(
		FILE * mesh_file,
		Eigen::PlainObjectBase< DerivedV > & V,
		Eigen::PlainObjectBase< DerivedT > & T,
		Eigen::PlainObjectBase< DerivedF > & F,
		Eigen::PlainObjectBase< DerivedF > & SF,
		Eigen::PlainObjectBase< DerivedC > & C,
		Eigen::PlainObjectBase< DerivedN > & N,
		Eigen::VectorXi & A)
{
	V.resize(0, Eigen::NoChange);
	T.resize(0, Eigen::NoChange);
	F.resize(0, Eigen::NoChange);
	SF.resize(0, Eigen::NoChange);
	C.resize(0, Eigen::NoChange);
	N.resize(0, Eigen::NoChange);
	A.resize(0, Eigen::NoChange);

	using namespace std;
#ifndef LINE_MAX
#  define LINE_MAX 2048
#endif
	char line[LINE_MAX];
	bool still_comments;

	// eat comments at beginning of file
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	char str[LINE_MAX];
	sscanf(line, " %s", str);
	// check that first word is MeshVersionFormatted
	if (0 != strcmp(str, "MeshVersionFormatted"))
	{
		fprintf(stderr,
			"Error: first word should be MeshVersionFormatted not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int one = -1;
	if (2 != sscanf(line, "%s %d", str, &one))
	{
		// 1 appears on next line?
		fscanf(mesh_file, " %d", &one);
	}
	if (one != 1)
	{
		fprintf(stderr, "Error: second word should be 1 not %d\n", one);
		fclose(mesh_file);
		return false;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that third word is Dimension
	if (0 != strcmp(str, "Dimension"))
	{
		fprintf(stderr, "Error: third word should be Dimension not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int three = -1;
	if (2 != sscanf(line, "%s %d", str, &three))
	{
		// 1 appears on next line?
		fscanf(mesh_file, " %d", &three);
	}
	if (three != 3)
	{
		fprintf(stderr, "Error: only Dimension 3 supported not %d\n", three);
		fclose(mesh_file);
		return false;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that fifth word is Vertices
	if (0 != strcmp(str, "Vertices"))
	{
		fprintf(stderr, "Error: fifth word should be Vertices not %s\n", str);
		fclose(mesh_file);
		return false;
	}

	//fgets(line,LINE_MAX,mesh_file);

	int number_of_vertices;
	if (1 != fscanf(mesh_file, " %d", &number_of_vertices) || number_of_vertices > 1000000000)
	{
		fprintf(stderr, "Error: expecting number of vertices less than 10^9...\n");
		fclose(mesh_file);
		return false;
	}
	// allocate space for vertices
	V.resize(number_of_vertices, 3);
	int extra;
	for (int i = 0; i<number_of_vertices; i++)
	{
		double x, y, z;
		if (4 != fscanf(mesh_file, " %lg %lg %lg %d", &x, &y, &z, &extra))
		{
			fprintf(stderr, "Error: expecting vertex position...\n");
			fclose(mesh_file);
			return false;
		}
		V(i, 0) = x;
		V(i, 1) = y;
		V(i, 2) = z;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that sixth word is Triangles
	if (0 != strcmp(str, "Triangles"))
	{
		fprintf(stderr, "Error: sixth word should be Triangles not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int number_of_triangles;
	if (1 != fscanf(mesh_file, " %d", &number_of_triangles))
	{
		fprintf(stderr, "Error: expecting number of triangles...\n");
		fclose(mesh_file);
		return false;
	}

	int number_of_surface_triangles = 0;
	SF.resize(number_of_triangles, 3);
	// triangle indices
	int tri[3];
	for (int i = 0; i<number_of_triangles; i++)
	{
		if (4 != fscanf(mesh_file, " %d %d %d %d", &tri[0], &tri[1], &tri[2], &extra))
		{
			printf("Error: expecting triangle indices...\n");
			return false;
		}
		if (extra == 1) { // 1 means boundary, 0 means non-boundary
			SF(number_of_surface_triangles, 0) = tri[0] - 1;
			SF(number_of_surface_triangles, 1) = tri[1] - 1;
			SF(number_of_surface_triangles, 2) = tri[2] - 1;
			number_of_surface_triangles++;
		}
	}
	SF.conservativeResize(number_of_surface_triangles, Eigen::NoChange);

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that sixth word is Triangles
	if (0 != strcmp(str, "Tetrahedra"))
	{
		fprintf(stderr, "Error: seventh word should be Tetrahedra not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int number_of_tetrahedra;
	if (1 != fscanf(mesh_file, " %d", &number_of_tetrahedra))
	{
		fprintf(stderr, "Error: expecting number of tetrahedra...\n");
		fclose(mesh_file);
		return false;
	}
	// allocate space for tetrahedra
	A.resize(number_of_tetrahedra); // region attribute
	T.resize(number_of_tetrahedra, 4);
	// tet indices
	int a, b, c, d;
	//[region]->[vertex index]
	// region == 0 => elastic tet, flesh tet
	// region < 0  => rigid tet, bone tet
	std::map<int, std::set<int>> rigid_vertex_map;
	std::unordered_set<int> rigid_vertex_set;
	for (int i = 0; i<number_of_tetrahedra; i++)
	{
		if (5 != fscanf(mesh_file, " %d %d %d %d %d", &a, &b, &c, &d, &extra))
		{
			fprintf(stderr, "Error: expecting tetrahedra indices...\n");
			fclose(mesh_file);
			return false;
		}
		T(i, 0) = a - 1;
		T(i, 1) = b - 1;
		T(i, 2) = c - 1;
		T(i, 3) = d - 1;
		A(i) = extra;

		if (extra < 0) {
			rigid_vertex_map[extra].insert(T(i, 0));
			rigid_vertex_map[extra].insert(T(i, 1));
			rigid_vertex_map[extra].insert(T(i, 2));
			rigid_vertex_map[extra].insert(T(i, 3));

			rigid_vertex_set.insert(T(i, 0));
			rigid_vertex_set.insert(T(i, 1));
			rigid_vertex_set.insert(T(i, 2));
			rigid_vertex_set.insert(T(i, 3));
		}
	}
	fclose(mesh_file);

	// rearrage the order of the vertices so that V = [Vf; Vt]
	std::unordered_map<int, int> RI;

	N.resize(rigid_vertex_map.size() + 1, Eigen::NoChange);
	int idx = 1;
	int nf = number_of_vertices;
	for (auto it = rigid_vertex_map.begin(); it != rigid_vertex_map.end(); it++)
	{
		N(idx) = it->second.size();
		nf -= N(idx);
		idx++;
	}
	N(0) = nf;

	Eigen::PlainObjectBase<DerivedV> rearranged_V;
	rearranged_V.resize(number_of_vertices, 3);
	int number_of_vertices_prev_set = 0;
	for (auto m_it = rigid_vertex_map.begin(); m_it != rigid_vertex_map.end(); m_it++)
	{
		if (m_it != rigid_vertex_map.begin()) {
			number_of_vertices_prev_set += std::prev(m_it)->second.size();
		}
		std::set<int> &one_rigid_set = m_it->second;
		for (auto s_it = one_rigid_set.begin(); s_it != one_rigid_set.end(); s_it++)
		{
			int cnt = number_of_vertices_prev_set + std::distance(one_rigid_set.begin(), s_it);
			rearranged_V.row(nf + cnt) = V.row(*s_it);
			RI[*s_it] = nf + cnt;
		}
	}

	int elastic_idx = 0;
	for (int i = 0; i < number_of_vertices; i++) {
		if (rigid_vertex_set.find(i) == rigid_vertex_set.end()) //elastic vertex
		{
			rearranged_V.row(elastic_idx) = V.row(i);
			RI[i] = elastic_idx;
			elastic_idx++;
		}
	}
	assert(elastic_idx == nf);
	V = rearranged_V;

	// update index of tetrahedron because V is rearranged
	for (int i = 0; i < number_of_tetrahedra; i++) {
		for (int j = 0; j < 4; j++) {
			T(i, j) = RI[T(i, j)];
		}
	}

	// update index of triangles because V is rearranged
	for (int i = 0; i < number_of_surface_triangles; i++) {
		for (int j = 0; j < 3; j++) {
			SF(i, j) = RI[SF(i, j)];
		}
	}

	// construct Triangle from Tetrahedron
	F.resize(4 * number_of_tetrahedra, 3);
	for (int i = 0; i < number_of_tetrahedra; ++i) {
		F(i * 4 + 0, 0) = T(i, 0); F(i * 4 + 0, 1) = T(i, 1); F(i * 4 + 0, 2) = T(i, 3);
		F(i * 4 + 1, 0) = T(i, 1); F(i * 4 + 1, 1) = T(i, 2); F(i * 4 + 1, 2) = T(i, 3);
		F(i * 4 + 2, 0) = T(i, 2); F(i * 4 + 2, 1) = T(i, 0); F(i * 4 + 2, 2) = T(i, 3);
		F(i * 4 + 3, 0) = T(i, 0); F(i * 4 + 3, 1) = T(i, 2); F(i * 4 + 3, 2) = T(i, 1);
	}

	// assign different colors to faces belong to rigid tets

	// http://colorbrewer2.org
	// 12 data classes
	// qualitative color scheme
	Eigen::MatrixX4d colorScheme;
	colorScheme.resize(12, 4);

	colorScheme <<
		166, 206, 227, 255
		, 31, 120, 180, 255
		, 178, 223, 138, 255
		, 51, 160, 44, 255
		, 251, 154, 153, 255
		, 227, 26, 28, 255
		, 253, 191, 111, 255
		, 255, 127, 0, 255
		, 202, 178, 214, 255
		, 106, 61, 154, 255
		, 255, 255, 153, 255
		, 177, 89, 40, 255;

	colorScheme /= 255.0f;
	C.resize(F.rows(), 4);
	for (int f_i = 0; f_i < F.rows(); ++f_i) {
		const int colorIdx = abs(A(f_i / 4)) % 12;
		Eigen::Vector4d color = colorScheme.row(colorIdx);
		C.row(f_i) = color;
	}
	return true;
}


// readMESH for self collision
template <
	typename DerivedV,
	typename DerivedT,
	typename DerivedF,
	typename DerivedC,
	typename DerivedN>
	IGL_INLINE bool igl::readMESH(
		const std::string mesh_file_name,
		Eigen::PlainObjectBase<DerivedV> & V,
		Eigen::PlainObjectBase<DerivedT> & T,
		Eigen::PlainObjectBase<DerivedF> & F,
		Eigen::PlainObjectBase<DerivedT> & ST,
		Eigen::PlainObjectBase<DerivedF> & SF,
		Eigen::PlainObjectBase<DerivedC> & C,
		Eigen::PlainObjectBase<DerivedN> & N,
		Eigen::VectorXi & A,
		Eigen::VectorXi & SV)
{
	using namespace std;
	FILE * mesh_file = fopen(mesh_file_name.c_str(), "r");
	if (NULL == mesh_file)
	{
		fprintf(stderr, "IOError: %s could not be opened...", mesh_file_name.c_str());
		return false;
	}
	return readMESH(mesh_file, V, T, F, ST, SF, C, N, A, SV);
}


template <
	typename DerivedV,
	typename DerivedT,
	typename DerivedF,
	typename DerivedC,
	typename DerivedN>
	IGL_INLINE bool igl::readMESH(
		FILE * mesh_file,
		Eigen::PlainObjectBase<DerivedV> & V,
		Eigen::PlainObjectBase<DerivedT> & T,
		Eigen::PlainObjectBase<DerivedF> & F,
		Eigen::PlainObjectBase<DerivedT> & ST,
		Eigen::PlainObjectBase<DerivedF> & SF,
		Eigen::PlainObjectBase<DerivedC> & C,
		Eigen::PlainObjectBase<DerivedN> & N,
		Eigen::VectorXi & A,
		Eigen::VectorXi & SV)
{
	V.resize(0, Eigen::NoChange);
	T.resize(0, Eigen::NoChange);
	F.resize(0, Eigen::NoChange);
	ST.resize(0, Eigen::NoChange);
	SF.resize(0, Eigen::NoChange);
	C.resize(0, Eigen::NoChange);
	N.resize(0, Eigen::NoChange);
	A.resize(0, Eigen::NoChange);
	SV.resize(0, Eigen::NoChange);

	using namespace std;
#ifndef LINE_MAX
#  define LINE_MAX 2048
#endif
	char line[LINE_MAX];
	bool still_comments;

	// eat comments at beginning of file
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	char str[LINE_MAX];
	sscanf(line, " %s", str);
	// check that first word is MeshVersionFormatted
	if (0 != strcmp(str, "MeshVersionFormatted"))
	{
		fprintf(stderr,
			"Error: first word should be MeshVersionFormatted not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int one = -1;
	if (2 != sscanf(line, "%s %d", str, &one))
	{
		// 1 appears on next line?
		fscanf(mesh_file, " %d", &one);
	}
	if (one != 1)
	{
		fprintf(stderr, "Error: second word should be 1 not %d\n", one);
		fclose(mesh_file);
		return false;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that third word is Dimension
	if (0 != strcmp(str, "Dimension"))
	{
		fprintf(stderr, "Error: third word should be Dimension not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int three = -1;
	if (2 != sscanf(line, "%s %d", str, &three))
	{
		// 1 appears on next line?
		fscanf(mesh_file, " %d", &three);
	}
	if (three != 3)
	{
		fprintf(stderr, "Error: only Dimension 3 supported not %d\n", three);
		fclose(mesh_file);
		return false;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that fifth word is Vertices
	if (0 != strcmp(str, "Vertices"))
	{
		fprintf(stderr, "Error: fifth word should be Vertices not %s\n", str);
		fclose(mesh_file);
		return false;
	}

	//fgets(line,LINE_MAX,mesh_file);

	int number_of_vertices;
	if (1 != fscanf(mesh_file, " %d", &number_of_vertices) || number_of_vertices > 1000000000)
	{
		fprintf(stderr, "Error: expecting number of vertices less than 10^9...\n");
		fclose(mesh_file);
		return false;
	}
	// allocate space for vertices
	V.resize(number_of_vertices, 3);
	int extra;
	for (int i = 0; i<number_of_vertices; i++)
	{
		double x, y, z;
		if (4 != fscanf(mesh_file, " %lg %lg %lg %d", &x, &y, &z, &extra))
		{
			fprintf(stderr, "Error: expecting vertex position...\n");
			fclose(mesh_file);
			return false;
		}
		V(i, 0) = x;
		V(i, 1) = y;
		V(i, 2) = z;
	}

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that sixth word is Triangles
	if (0 != strcmp(str, "Triangles"))
	{
		fprintf(stderr, "Error: sixth word should be Triangles not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int number_of_triangles;
	if (1 != fscanf(mesh_file, " %d", &number_of_triangles))
	{
		fprintf(stderr, "Error: expecting number of triangles...\n");
		fclose(mesh_file);
		return false;
	}

	// triangle indices
	int tri[3];
	std::set<int> SVI; // surface vertex index
	int number_of_surface_triangles = 0;
	SF.resize(number_of_triangles, 3);
	for (int i = 0; i<number_of_triangles; i++)
	{
		if (4 != fscanf(mesh_file, " %d %d %d %d", &tri[0], &tri[1], &tri[2], &extra))
		{
			printf("Error: expecting triangle indices...\n");
			return false;
		}

		if (extra == 1) { // 1 means boundary, 0 means non-boundary
			SVI.insert(tri[0] - 1);
			SVI.insert(tri[1] - 1);
			SVI.insert(tri[2] - 1);
			SF(number_of_surface_triangles, 0) = tri[0] - 1;
			SF(number_of_surface_triangles, 1) = tri[1] - 1;
			SF(number_of_surface_triangles, 2) = tri[2] - 1;
			number_of_surface_triangles++;
		}
	}
	SF.conservativeResize(number_of_surface_triangles, 3);

	std::cout << "SF = " << SF << std::endl;

	// eat comments
	still_comments = true;
	while (still_comments)
	{
		fgets(line, LINE_MAX, mesh_file);
		still_comments = (line[0] == '#' || line[0] == '\n');
	}

	sscanf(line, " %s", str);
	// check that sixth word is Triangles
	if (0 != strcmp(str, "Tetrahedra"))
	{
		fprintf(stderr, "Error: seventh word should be Tetrahedra not %s\n", str);
		fclose(mesh_file);
		return false;
	}
	int number_of_tetrahedra;
	if (1 != fscanf(mesh_file, " %d", &number_of_tetrahedra))
	{
		fprintf(stderr, "Error: expecting number of tetrahedra...\n");
		fclose(mesh_file);
		return false;
	}
	// allocate space for tetrahedra
	A.resize(number_of_tetrahedra); // region attribute
	T.resize(number_of_tetrahedra, 4);
	// tet indices
	int a, b, c, d;
	//[region]->[vertex index]
	// region == 0 => elastic tet, flesh tet
	// (region < 0 && region % 10 == 0) => rigid tet, bone tet (e.g. -10, -20, -30)
	// (region > 0 && region % 10 == 0) => activation tet, muscle tet (e.g. 10, 20, 30)
	std::map<int, std::set<int>> rigid_vertex_map;
	std::map<int, std::set<int>> active_vertex_map;
	std::unordered_set<int> rigid_vertex_set;
	std::unordered_set<int> active_vertex_set;
	int number_of_surface_tetrahedra = 0;
	for (int i = 0; i<number_of_tetrahedra; i++)
	{
		if (5 != fscanf(mesh_file, " %d %d %d %d %d", &a, &b, &c, &d, &extra))
		{
			fprintf(stderr, "Error: expecting tetrahedra indices...\n");
			fclose(mesh_file);
			return false;
		}
		T(i, 0) = a - 1;
		T(i, 1) = b - 1;
		T(i, 2) = c - 1;
		T(i, 3) = d - 1;
		A(i) = extra;

		if (extra < 0) {
			rigid_vertex_map[extra].insert(T(i, 0));
			rigid_vertex_map[extra].insert(T(i, 1));
			rigid_vertex_map[extra].insert(T(i, 2));
			rigid_vertex_map[extra].insert(T(i, 3));

			rigid_vertex_set.insert(T(i, 0));
			rigid_vertex_set.insert(T(i, 1));
			rigid_vertex_set.insert(T(i, 2));
			rigid_vertex_set.insert(T(i, 3));
		}
		else if (extra % 10 == 0) {
			active_vertex_map[extra].insert(T(i, 0));
			active_vertex_map[extra].insert(T(i, 1));
			active_vertex_map[extra].insert(T(i, 2));
			active_vertex_map[extra].insert(T(i, 3));

			active_vertex_set.insert(T(i, 0));
			active_vertex_set.insert(T(i, 1));
			active_vertex_set.insert(T(i, 2));
			active_vertex_set.insert(T(i, 3));
		}
	}
	fclose(mesh_file);

	// rearrage the order of the vertices so that V = [Vf; Vt]
	std::unordered_map<int, int> RI;

	N.resize(rigid_vertex_map.size() + 1, Eigen::NoChange);
	int idx = 1;
	int nf = number_of_vertices;
	for (auto it = rigid_vertex_map.begin(); it != rigid_vertex_map.end(); it++)
	{
		N(idx) = it->second.size();
		nf -= N(idx);
		idx++;
	}
	N(0) = nf;

	Eigen::PlainObjectBase<DerivedV> rearranged_V;
	rearranged_V.resize(number_of_vertices, 3);
	int number_of_vertices_prev_set = 0;
	for (auto m_it = rigid_vertex_map.begin(); m_it != rigid_vertex_map.end(); m_it++)
	{
		if (m_it != rigid_vertex_map.begin()) {
			number_of_vertices_prev_set += std::prev(m_it)->second.size();
		}
		std::set<int> &one_rigid_set = m_it->second;
		for (auto s_it = one_rigid_set.begin(); s_it != one_rigid_set.end(); s_it++)
		{
			int cnt = number_of_vertices_prev_set + std::distance(one_rigid_set.begin(), s_it);
			rearranged_V.row(nf + cnt) = V.row(*s_it);
			RI[*s_it] = nf + cnt;
		}
	}

	int elastic_idx = 0;
	for (int i = 0; i < number_of_vertices; i++) {
		if (rigid_vertex_set.find(i) == rigid_vertex_set.end()) //elastic vertex
		{
			rearranged_V.row(elastic_idx) = V.row(i);
			RI[i] = elastic_idx;
			elastic_idx++;
		}
	}
	assert(elastic_idx == nf);
	V = rearranged_V;

	// update index of tetrahedron because V is rearranged
	for (int i = 0; i < number_of_tetrahedra; i++) {
		if (SVI.find(T(i, 0)) != SVI.end() ||
			SVI.find(T(i, 1)) != SVI.end() ||
			SVI.find(T(i, 2)) != SVI.end() ||
			SVI.find(T(i, 3)) != SVI.end()) { // tetrahedron on surface
			ST.conservativeResize(number_of_surface_tetrahedra + 1, 4);
			ST(number_of_surface_tetrahedra, 0) = RI[T(i, 0)];
			ST(number_of_surface_tetrahedra, 1) = RI[T(i, 1)];
			ST(number_of_surface_tetrahedra, 2) = RI[T(i, 2)];
			ST(number_of_surface_tetrahedra, 3) = RI[T(i, 3)];
			number_of_surface_tetrahedra++;
		}

		for (int j = 0; j < 4; j++) {
			T(i, j) = RI[T(i, j)];
		}
	}

	SV.resize(SVI.size());
	int number_of_surface_vertices = 0;
	for (auto it = SVI.begin(); it != SVI.end(); it++) {
		SV(number_of_surface_vertices) = RI[*it];
		number_of_surface_vertices++;
	}

	// construct Triangle from Tetrahedron
	F.resize(4 * number_of_tetrahedra, 3);
	for (int i = 0; i < number_of_tetrahedra; ++i) {
		F(i * 4 + 0, 0) = T(i, 0); F(i * 4 + 0, 1) = T(i, 1); F(i * 4 + 0, 2) = T(i, 3);
		F(i * 4 + 1, 0) = T(i, 1); F(i * 4 + 1, 1) = T(i, 2); F(i * 4 + 1, 2) = T(i, 3);
		F(i * 4 + 2, 0) = T(i, 2); F(i * 4 + 2, 1) = T(i, 0); F(i * 4 + 2, 2) = T(i, 3);
		F(i * 4 + 3, 0) = T(i, 0); F(i * 4 + 3, 1) = T(i, 2); F(i * 4 + 3, 2) = T(i, 1);
	}

	// assign different colors to faces belong to rigid tets

	// http://colorbrewer2.org
	// 12 data classes
	// qualitative color scheme
	Eigen::MatrixX4d colorScheme;
	colorScheme.resize(12, 4);

	colorScheme <<
		166, 206, 227, 255
		, 31, 120, 180, 255
		, 178, 223, 138, 255
		, 51, 160, 44, 255
		, 251, 154, 153, 255
		, 227, 26, 28, 255
		, 253, 191, 111, 255
		, 255, 127, 0, 255
		, 202, 178, 214, 255
		, 106, 61, 154, 255
		, 255, 255, 153, 255
		, 177, 89, 40, 255;

	colorScheme /= 255.0f;
	C.resize(F.rows(), 4);
	for (int f_i = 0; f_i < F.rows(); ++f_i) {
		const int colorIdx = abs(A(f_i / 4)) % 12;
		Eigen::Vector4d color = colorScheme.row(colorIdx);
		C.row(f_i) = color;
	}
	return true;
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template bool igl::readMESH<Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(FILE*, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
// generated by autoexplicit.sh
template bool igl::readMESH<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(FILE*, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> >&);
// generated by autoexplicit.sh
template bool igl::readMESH<Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template bool igl::readMESH<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<unsigned int, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, 3, 1, -1, 3> >&);
// generated by autoexplicit.sh
template bool igl::readMESH<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
// generated by autoexplicit.sh
template bool igl::readMESH<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template bool igl::readMESH<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template bool igl::readMESH<double, int>(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >&);
template bool igl::readMESH< Eigen::Matrix<double, -1, -1, 0, -1, -1>,  Eigen::Matrix<int, -1, -1, 0, -1, -1>,  Eigen::Matrix<int, -1, -1, 0, -1, -1>,  Eigen::Matrix<double, -1, -1, 0, -1, -1>,  Eigen::Matrix<int, -1, 1, 0, -1, 1> >( std::basic_string<char, struct std::char_traits<char>,  std::allocator<char> >,  Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > &,  Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > &,  Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > &,  Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > &,  Eigen::PlainObjectBase< Eigen::Matrix<int, -1, 1, 0, -1, 1> > &,  Eigen::Matrix<int, -1, 1, 0, -1, 1> &);
template bool igl::readMESH< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(std::basic_string<char, struct std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, 1, 0, -1, 1> > &, Eigen::Matrix<int, -1, 1, 0, -1, 1> &);
#endif