#include "self_collision.h"
#include "spatial_hash.h"
#include <Eigen/Geometry>
#include <igl/per_face_normals.h>
#include <igl/signed_distance.h>
#include <igl/volume.h>
#include <iostream>
#include <algorithm>


extern double s_grid;


template <typename Scalar>
IGL_INLINE bool igl::point_in_tet(
	const Eigen::RowVector3d & vtx,
	const Eigen::Matrix<Scalar, 4, 3> & tet
)
{
	using namespace Eigen;

	double x = vtx.x();
	double y = vtx.y();
	double z = vtx.z();

	double x1 = tet(0, 0);
	double y1 = tet(0, 1);
	double z1 = tet(0, 2);

	double x2 = tet(1, 0);
	double y2 = tet(1, 1);
	double z2 = tet(1, 2);

	double x3 = tet(2, 0);
	double y3 = tet(2, 1);
	double z3 = tet(2, 2);

	double x4 = tet(3, 0);
	double y4 = tet(3, 1);
	double z4 = tet(3, 2);


	//http://steve.hollasch.net/cgindex/geometry/ptintet.html
	Matrix4d D0, D1, D2, D3, D4;
	D0 << x1, y1, z1, 1,
		x2, y2, z2, 1,
		x3, y3, z3, 1,
		x4, y4, z4, 1;

	D1 << x, y, z, 1,
		x2, y2, z2, 1,
		x3, y3, z3, 1,
		x4, y4, z4, 1;

	D2 << x1, y1, z1, 1,
		x, y, z, 1,
		x3, y3, z3, 1,
		x4, y4, z4, 1;

	D3 << x1, y1, z1, 1,
		x2, y2, z2, 1,
		x, y, z, 1,
		x4, y4, z4, 1;

	D4 << x1, y1, z1, 1,
		x2, y2, z2, 1,
		x3, y3, z3, 1,
		x, y, z, 1;

	double eps = 1e-8;
	double det[5] = { D0.determinant(), D1.determinant(), D2.determinant(), D3.determinant(), D4.determinant() };

	if (det[0] > eps) { // inverted tet
		return false;
	}
	if (det[0] < -eps && det[1] < -eps && det[2] < -eps && det[3] < -eps && det[4] < -eps) {
		return true;
	}
	return false;
}


template <typename Scalar>
IGL_INLINE bool igl::edge_surface_intersection(
	const Eigen::Matrix<Scalar, 2, 3> & edge,
	const Eigen::Matrix<Scalar, 3, 3> & tri,
	Scalar &t)
{
	using namespace Eigen;

	Vector3d p1    = tri.row(0).transpose();
	Vector3d p2    = tri.row(1).transpose();
	Vector3d p3    = tri.row(2).transpose();
	Vector3d p     = edge.row(0).transpose();
	Vector3d p_end = edge.row(1).transpose();

	Vector3d e  = p - p1;
	Vector3d e1 = p2 - p1;
	Vector3d e2 = p3 - p1;
	Vector3d d = p_end - p;

	Vector3d result = Vector3d(
		d.cross(e2).dot(e),
		e.cross(e1).dot(d),
		e.cross(e1).dot(e2));

	result /= d.cross(e2).dot(e1);

	Scalar b1 = result(0);
	Scalar b2 = result(1);
	t  = result(2);

	double eps = 1e-8;
	if (b1 >= 0 && b2 >= 0 && b1 + b2 <= 1 && t > eps && t < 1 - eps) {
		return true;
	}
	return false;
}



template <typename DerivedT>
	IGL_INLINE void igl::find_neighbors(
		const Eigen::PlainObjectBase< DerivedT > & T,
		std::map<int, std::set<int>> & neighbors)
{
	for (int tet_i = 0; tet_i < T.rows(); tet_i++) {
		for (int i = 0; i < 4; i++) {
			neighbors[T(tet_i, i)].insert(T(tet_i, (i + 1) % 4));
			neighbors[T(tet_i, i)].insert(T(tet_i, (i + 2) % 4));
			neighbors[T(tet_i, i)].insert(T(tet_i, (i + 3) % 4));
		}
	}
}


template <typename DerivedV, typename DerivedT>
IGL_INLINE double igl::compute_grid_size(
	const Eigen::PlainObjectBase< DerivedV > & V,
	const Eigen::PlainObjectBase< DerivedT > & T
)
{
	double grid_size = 0.; // grid size is the average size of a tetrahedron
	for (int tet_i = 0; tet_i < T.rows(); tet_i++) {
		for (int i = 0; i < 4; i++) {
			grid_size += (V.row(T(tet_i, 0)) - V.row(T(tet_i, 3))).norm();
			grid_size += (V.row(T(tet_i, 1)) - V.row(T(tet_i, 3))).norm();
			grid_size += (V.row(T(tet_i, 2)) - V.row(T(tet_i, 3))).norm();
		}
	}
	grid_size /= (12 * T.rows());
	std::cout << "grid_size = " << grid_size << std::endl;
	return grid_size;
}

#ifdef IGL_STATIC_LIBRARY
template bool igl::self_collision< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::Matrix<int, -1, 1, 0, -1, 1> const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::vector< Eigen::Triplet<double, int>, std::allocator< Eigen::Triplet<double, int> > > &);
template bool igl::self_collision< Eigen::Matrix<double, -1, -1, 0, -1, -1>,  Eigen::Matrix<int, -1, -1, 0, -1, -1>,  Eigen::Matrix<int, -1, -1, 0, -1, -1> >( Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,  Eigen::Matrix<int, -1, 1, 0, -1, 1> const &,  Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,  Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,  std::set<int, struct std::less<int>,  std::allocator<int> > const &,  Eigen::Matrix<int, -1, 1, 0, -1, 1> &,  Eigen::Matrix<double, -1, 3, 0, -1, 3> &);
template bool igl::point_collisions< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::set<int, struct std::less<int>, std::allocator<int> > const &, std::set<int, struct std::less<int>, std::allocator<int> > &);
template void igl::find_neighbors< Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::map<int, std::set<int, struct std::less<int>, std::allocator<int> >, struct std::less<int>, std::allocator<struct std::pair<int const, std::set<int, struct std::less<int>, std::allocator<int> > > > > &);
template bool igl::self_collision< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::map<int, std::set<int, struct std::less<int>, std::allocator<int> >, struct std::less<int>, std::allocator<struct std::pair<int const, std::set<int, struct std::less<int>, std::allocator<int> > > > > const &, std::set<int, struct std::less<int>, std::allocator<int> > const &, struct igl::RBCData &);
template double igl::compute_grid_size< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &);
template bool igl::point_in_tet<double>(Eigen::Matrix<double, 1, 3, 1, 1, 3> const &, Eigen::Matrix<double, 4, 3, 0, 4, 3> const &);
template bool igl::edge_surface_intersection<double>(Eigen::Matrix<double, 2, 3, 0, 2, 3> const &, Eigen::Matrix<double, 3, 3, 0, 3, 3> const &, double &);


#endif