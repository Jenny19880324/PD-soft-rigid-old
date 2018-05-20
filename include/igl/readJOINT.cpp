// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "readJOINT.h"
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <iostream>

template <typename Index, typename DerivedV>
IGL_INLINE bool igl::readJOINT(
const std::string joint_file_name,
std::vector<std::vector<Index>> & I,
Eigen::PlainObjectBase<DerivedV> & V)
{
  using namespace std;
  FILE * joint_file = fopen(joint_file_name.c_str(),"r");
  if(NULL==joint_file)
  {
    fprintf(stderr,"IOError: %s could not be opened...",joint_file_name.c_str());
    return false;
  }
  return igl::readJOINT(joint_file,I, V);
}


template <typename Index, typename DerivedV>
IGL_INLINE bool igl::readJOINT(
FILE * joint_file,
std::vector<std::vector<Index>> & I,
Eigen::PlainObjectBase<DerivedV> & V)
{
  using namespace std;
#ifndef LINE_MAX
#  define LINE_MAX 2048
#endif
  char line[LINE_MAX];
  bool still_comments;
  I.clear();

  // eat comments at beginning of file
  still_comments= true;
  while(still_comments)
  {
    if(fgets(line,LINE_MAX,joint_file) == NULL)
    {
      fprintf(stderr, "Error: couldn't find start of .joint file");
      fclose(joint_file);
      return false;
    }
    still_comments = (line[0] == '#' || line[0] == '\n');
  }
  char str[LINE_MAX];

  sscanf(line," %s",str);
  // check that the first word is Joints
  if(0!=strcmp(str,"Joints"))
  {
    fprintf(stderr,"Error: first word should be Joints not %s\n",str);
    fclose(joint_file);
    return false;
  }

  int number_of_joints;
  if(1 != fscanf(joint_file," %d",&number_of_joints) || number_of_joints > 1000000000)
  {
    fprintf(stderr,"Error: expecting number of joints less than 10^9...\n");
    fclose(joint_file);
    return false;
  }
  // allocate space for joints
  I.resize(number_of_joints);
  V.resize(number_of_joints, 3);
  for(int i = 0;i<number_of_joints;i++)
  {
	  int number_of_rigid_bodies_involved;
	  if (1 != fscanf(joint_file, "%d", &number_of_rigid_bodies_involved) || number_of_rigid_bodies_involved > 1000000000)
	  {
		  fprintf(stderr, "Error: expecting number of rigid body involved less than 10^9...\n");
		  fclose(joint_file);
		  return false;
	  }
	  std::vector<Index> idx_of_one_joint;
	  for (int j = 0; j < number_of_rigid_bodies_involved; j++) 
	  {
		  int idx;
		  if (1 != fscanf(joint_file, "%d", &idx)) {
			  fprintf(stderr, "Error: expecting index of rigid bodies, got %d\n", idx);
			  fclose(joint_file);
			  return false;
		  }
		  idx_of_one_joint.push_back(idx);
	  }
	  I[i] = idx_of_one_joint;
	  
    double x,y,z;
    if(3 != fscanf(joint_file," %lg %lg %lg",&x,&y,&z))
    {
      fprintf(stderr,"Error: expecting vertex position...\n");
      fclose(joint_file);
      return false;
    }
	V(i, 0) = x;
	V(i, 1) = y;
	V(i, 2) = z;
  }
  fclose(joint_file);
  return true;
}


#ifdef IGL_STATIC_LIBRARY
template bool igl::readJOINT<int, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, struct std::char_traits<char>, std::allocator<char> >, std::vector<class std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);
#endif