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

template <typename DerivedV, typename DerivedT, typename DerivedF, typename DerivedC>
IGL_INLINE bool igl::readMESH(
	const std::string mesh_file_name,
	Eigen::PlainObjectBase<DerivedV>& V,
	Eigen::PlainObjectBase<DerivedV>& Vf,
	Eigen::PlainObjectBase<DerivedV>& Vb,
	Eigen::PlainObjectBase<DerivedT>& T,
	Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedC>& C)
{
	using namespace std;
	FILE * mesh_file = fopen(mesh_file_name.c_str(), "r");
	if (NULL == mesh_file)
	{
		fprintf(stderr, "IOError: %s could not be opened...", mesh_file_name.c_str());
		return false;
	}
	return readMESH(mesh_file, V, Vf, Vb, T, F, C);
}


template <typename DerivedV, typename DerivedT, typename DerivedF, typename DerivedC>
IGL_INLINE bool igl::readMESH(
	FILE * mesh_file,
	Eigen::PlainObjectBase<DerivedV>& V,
	Eigen::PlainObjectBase<DerivedV>& Vf,
	Eigen::PlainObjectBase<DerivedV>& Vb,
	Eigen::PlainObjectBase<DerivedT>& T,
	Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedC>& C)
{
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
	T.resize(number_of_tetrahedra, 4);
	// tet indices
	int a, b, c, d;
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

		if (extra == 0) {
			rigid_vertex_set.insert(T(i, 0));
			rigid_vertex_set.insert(T(i, 1));
			rigid_vertex_set.insert(T(i, 2));
			rigid_vertex_set.insert(T(i, 3));
		}
	}
	fclose(mesh_file);

	// rearrage the order of the vertices so that V = [Vf; Vt]
	std::unordered_map<int, int> RI;
	int nb = rigid_vertex_set.size();
	int nf = number_of_vertices - nb;
	for (int i = 0; i < number_of_vertices; i++) {
		if (rigid_vertex_set.find(i) == rigid_vertex_set.end()) { // flesh vertices
			RI[i] = Vf.rows();
			Vf.conservativeResize(Vf.rows() + 1, 3);
			Vf.row(RI[i]) = V.row(i);
		}
		else { // bone vertices
			RI[i] = nf + Vb.rows();
			Vb.conservativeResize(Vb.rows() + 1, 3);
			Vb.row(RI[i] - nf) = V.row(i);
		}
	}
	
	rigid_vertex_set.clear();
	for (int i = 0; i < nb; i++) {
		rigid_vertex_set.insert(nf + i);
	}
	assert(Vf.rows() == nf);
	assert(Vb.rows() == nb);
	V.setZero();
	V << Vf,
		Vb;

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
		bool is_rigid_face = rigid_vertex_set.find(F(f_i, 0)) != rigid_vertex_set.end() &&
			rigid_vertex_set.find(F(f_i, 1)) != rigid_vertex_set.end() &&
			rigid_vertex_set.find(F(f_i, 2)) != rigid_vertex_set.end();

		if (is_rigid_face) {
			const size_t colorIdx = 5;
			Eigen::Vector4d color = colorScheme.row(colorIdx);
			C.row(f_i) = color;
		}
		else {
			const size_t colorIdx = 1;
			Eigen::Vector4d color = colorScheme.row(colorIdx);
			C.row(f_i) = color;
		}
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
template bool igl::readMESH<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, struct std::char_traits<char>, class std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);
#endif