#include "self_collision.h"
#include "spatial_hash.h"
#include <Eigen/Geometry>
#include <igl/per_face_normals.h>
#include <igl/signed_distance.h>
#include <igl/volume.h>
#include <iostream>
#include <algorithm>


extern double s_grid;

template<typename DerivedV, typename DerivedT, typename Scalar>
IGL_INLINE bool igl::self_collision(
	const Eigen::PlainObjectBase <DerivedV> & V,
	const Eigen::VectorXi & SV,
	const Eigen::PlainObjectBase < DerivedT > & ST,
	std::vector<Eigen::Triplet<Scalar>> & L_triplets)
{
	using namespace std;
	using namespace Eigen;
	typedef Matrix<typename DerivedV::Scalar, 1, 3> RowVectorV3;

	for (int i = 0; i < SV.rows(); i++) {
		Scalar x = V(SV(i), 0);
		Scalar y = V(SV(i), 1);
		Scalar z = V(SV(i), 2);
		for (int tet_i = 0; tet_i < ST.rows(); tet_i++) {
			int t1 = ST(tet_i, 0);
			int t2 = ST(tet_i, 1);
			int t3 = ST(tet_i, 2);
			int t4 = ST(tet_i, 3);

			Scalar x1 = V(t1, 0);
			Scalar y1 = V(t1, 1);
			Scalar z1 = V(t1, 2);

			Scalar x2 = V(t2, 0);
			Scalar y2 = V(t2, 1);
			Scalar z2 = V(t2, 2);

			Scalar x3 = V(t3, 0);
			Scalar y3 = V(t3, 1);
			Scalar z3 = V(t3, 2);

			Scalar x4 = V(t4, 0);
			Scalar y4 = V(t4, 1);
			Scalar z4 = V(t4, 2);

			
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


			
			if (D0.determinant() > 0. && 
				D1.determinant() > 0. && 
				D2.determinant() > 0. && 
				D3.determinant() > 0. && 
				D4.determinant() > 0.){
				 std::cout << "vertex " << SV(i) << " collides with tet of vert " << t1 << ", " << t2 << ", " << t3 << ", " << t4 << std::endl;
				 std::cout << "v = " << x << ", " << y << ", " << z << std::endl;
				 std::cout << "D0 = " << D0 << std::endl;
				 std::cout << t1 << " " << t2 << " " << t3 << " " << t4 << std::endl;
			}
		}
	}
	return false;
}


template<
	typename DerivedV, 
	typename DerivedT,
	typename DerivedF>
IGL_INLINE bool igl::self_collision(
	const Eigen::PlainObjectBase < DerivedV > & V,
	const Eigen::VectorXi & SV,
	const Eigen::PlainObjectBase < DerivedT > & ST,
	const Eigen::PlainObjectBase < DerivedF > & SF,
	const std::set<int> &b_set,
	Eigen::VectorXi & b,
	Eigen::MatrixX3d & bc)
{
	//using namespace std;
	//using namespace Eigen;
	//typedef Matrix<typename DerivedV::Scalar, 1, 3> RowVectorV3;

	//b.resize(0);
	//bc.resize(0, Eigen::NoChange);

	//for (int i = 0; i < SV.rows(); i++) {
	//	if (b_set.find(SV(i)) != b_set.end()) {
	//		continue;
	//	}
	//	double x = V(SV(i), 0);
	//	double y = V(SV(i), 1);
	//	double z = V(SV(i), 2);
	//	for (int tet_i = 0; tet_i < ST.rows(); tet_i++) {
	//		int t1 = ST(tet_i, 0);
	//		int t2 = ST(tet_i, 1);
	//		int t3 = ST(tet_i, 2);
	//		int t4 = ST(tet_i, 3);

	//		double x1 = V(t1, 0);
	//		double y1 = V(t1, 1);
	//		double z1 = V(t1, 2);

	//		double x2 = V(t2, 0);
	//		double y2 = V(t2, 1);
	//		double z2 = V(t2, 2);

	//		double x3 = V(t3, 0);
	//		double y3 = V(t3, 1);
	//		double z3 = V(t3, 2);

	//		double x4 = V(t4, 0);
	//		double y4 = V(t4, 1);
	//		double z4 = V(t4, 2);


	//		//http://steve.hollasch.net/cgindex/geometry/ptintet.html
	//		Matrix4d D0, D1, D2, D3, D4;
	//		D0 << x1, y1, z1, 1,
	//			x2, y2, z2, 1,
	//			x3, y3, z3, 1,
	//			x4, y4, z4, 1;

	//		D1 << x, y, z, 1,
	//			x2, y2, z2, 1,
	//			x3, y3, z3, 1,
	//			x4, y4, z4, 1;

	//		D2 << x1, y1, z1, 1,
	//			x, y, z, 1,
	//			x3, y3, z3, 1,
	//			x4, y4, z4, 1;

	//		D3 << x1, y1, z1, 1,
	//			x2, y2, z2, 1,
	//			x, y, z, 1,
	//			x4, y4, z4, 1;

	//		D4 << x1, y1, z1, 1,
	//			x2, y2, z2, 1,
	//			x3, y3, z3, 1,
	//			x, y, z, 1;

	//		double eps = 1e-8;
	//		double det[5] = { D0.determinant(), D1.determinant(), D2.determinant(), D3.determinant(), D4.determinant() };

	//		if (det[0] > eps) {
	//			continue;
	//		}
	//		if (det[0] < -eps && det[1] < -eps && det[2] < -eps && det[3] < -eps && det[4] < -eps) {
	//			std::cout << "vertex " << SV(i) << " collides with tet of vert " << t1 << ", " << t2 << ", " << t3 << ", " << t4 << std::endl;
	//			
	//			// project a point to the triangle
	//			Vector3d p = Vector3d(x, y, z);
	//			Vector3d u = min_p2 - min_p1;
	//			Vector3d v = min_p3 - min_p1;
	//			Vector3d n = u.cross(v);
	//			Vector3d w = p - min_p1;
	//			Vector3d p_projected = p;

	//			double gamma = (u.cross(w)).dot(n) / (n.squaredNorm());
	//			double beta = (w.cross(v)).dot(n) / (n.squaredNorm());
	//			double alpha = 1. - gamma - beta;

	//			if (alpha >= 0 && alpha <= 1 &&
	//				beta >= 0 && beta <= 1 &&
	//				gamma >= 0 && gamma <= 1) {
	//				std::cout << "The point p lies inside T" << std::endl;
	//				p_projected = alpha * min_p1 + beta * min_p2 + gamma * min_p3;
	//				b.conservativeResize(b.rows() + 1);
	//				b(b.rows() - 1) = SV(i);
	//				bc.conservativeResize(b.rows(), Eigen::NoChange);
	//				bc.row(b.rows() - 1) = p_projected.transpose();

	//				//std::cout << "min_h = " << min_h << std::endl;
	//				//std::cout << "b = " << i << std::endl;
	//				//std::cout << "bc = " << p_projected << std::endl;
	//				//std::cout << "p = " << p << std::endl;
	//			}
	//		}
	//	}
	//}
	return false;
}


template<
	typename DerivedV,
	typename DerivedT,
	typename DerivedF>
	IGL_INLINE bool igl::self_collision(
		const Eigen::PlainObjectBase< DerivedV > & V,
		const Eigen::PlainObjectBase< DerivedT > & T,
		const Eigen::PlainObjectBase< DerivedF > & SF,
		const std::map<int, std::set<int>> & neighbors,
		const std::set<int> & b_set,
		RBCData & data)
{
	std::set<int> colliding_points;
	std::map<int, std::set<int>> border_points;
	std::set<igl::Edge> intersecting_edges;
	std::map<int, Eigen::Vector3d> border_penetration_direction;
	std::map<int, double> border_penetration_depth;
	std::map<int, Eigen::Vector3d> penetration_direction;
	std::map<int, double> penetration_depth;
	std::vector<Eigen::Vector3d> intersection_points;
	std::vector<Eigen::Vector3d> surface_normals;

	// debug
	std::set<int> non_collide_set = b_set;
	for (int i = 0; i < T.rows(); i++) {
		int t0 = T(i, 0);
		int t1 = T(i, 1);
		int t2 = T(i, 2);
		int t3 = T(i, 3);

		Eigen::Matrix<double, 4, 4> tet;
		tet << V.row(t0), 1,
			V.row(t1), 1,
			V.row(t2), 1,
			V.row(t3), 1;

		if (tet.determinant() > 1e-16) {
			non_collide_set.insert(t0);
			non_collide_set.insert(t1);
			non_collide_set.insert(t2);
			non_collide_set.insert(t3);
			//std::cout << "t0 = " << t0 << " t1 = " << t1 << " t2 = " << t2 << " t3 = " << t3 << std::endl;
		}
	}
	// debug
	// step 1
  	if (igl::point_collisions(V, T, non_collide_set, colliding_points)) {
		std::cout << "point collision" << std::endl;
	}

	// step 2
	if (igl::edge_intersections(V, T, SF, neighbors, colliding_points, non_collide_set, intersection_points, surface_normals, border_points)) {
		std::cout << "edge intersection" << std::endl;
	}

	// step 3
	if (igl::penetration_depth_and_direction(V, intersection_points, surface_normals, neighbors, border_points, border_penetration_direction, border_penetration_depth)) {
		std::cout << "penetration_depth_and_direction" << std::endl;
	}

	// step 4
	if (igl::propagation(V, border_penetration_direction, border_penetration_depth, colliding_points, border_points, neighbors, penetration_direction, penetration_depth)) {
		std::cout << "propagation" << std::endl;
	}

	//data.f_ext.resize(V.rows(), 3);
	//data.f_ext.setZero();

		std::cout << "data.self_collision_weight = " << data.self_collision_weight << std::endl;
	for (int v_i = 0; v_i < V.rows(); v_i++) {
		if (penetration_depth.find(v_i) != penetration_depth.end()){
			data.f_ext.row(v_i) += data.self_collision_weight * penetration_depth.at(v_i) * (penetration_direction.at(v_i).transpose());
		}
	}
	//data.f_ext += data.M * (Eigen::RowVector3d(0., (double)data.g, 0.).replicate(V.rows(), 1));
 	return !penetration_depth.empty();
}


template<typename DerivedV, typename DerivedT>
IGL_INLINE bool igl::point_collisions(
	const Eigen::PlainObjectBase< DerivedV > & V,
	const Eigen::PlainObjectBase< DerivedT > & T,
	const std::set<int> & b_set,
	std::set<int> & colliding_points
)
{
	// spatial hashing
	int number_of_buckets = T.rows();
	igl::SpatialHash<DerivedV> spatial_hash(s_grid, number_of_buckets);
	for (int v_i = 0; v_i < V.rows(); v_i++) {
		spatial_hash.insert(v_i);
	}
	for (int tet_i = 0; tet_i < T.rows(); tet_i++) {
		spatial_hash.insert(Tetrahedron(T(tet_i, 0), T(tet_i, 1), T(tet_i, 2), T(tet_i, 3), tet_i));
	}

	return spatial_hash.point_tet_collision(b_set, colliding_points);
}


template <
	typename DerivedV,
	typename DerivedT,
	typename DerivedF>
	IGL_INLINE bool igl::edge_intersections(
		const Eigen::PlainObjectBase< DerivedV > & V,
		const Eigen::PlainObjectBase< DerivedT > & T,
		const Eigen::PlainObjectBase< DerivedF > & SF,
		const std::map<int, std::set<int>> & neighbors,
		const std::set<int> & colliding_points,
		const std::set<int> & b_set,
		std::vector<Eigen::Vector3d> & intersection_points,
		std::vector<Eigen::Vector3d> & surface_normals,
		std::map<int, std::set<int>> & border_points)
{
	bool is_intersected = false;
	using namespace Eigen;
	std::set<Edge> intersecting_edges;
	assert(neighbors.size() == V.rows());
	int e_idx = 0;
	for (auto it = colliding_points.begin(); it != colliding_points.end(); it++) {
		if (b_set.find(*it) != b_set.end()) {
			continue;
		}
		for (auto n_it = neighbors.at(*it).begin(); n_it != neighbors.at(*it).end(); n_it++) {
			if (b_set.find(*n_it) != b_set.end()) {
				continue;
			}
			if (colliding_points.find(*n_it) == colliding_points.end()) {
				intersecting_edges.insert(Edge(*it, *n_it, e_idx));
				e_idx++;
			}
		}
	}

	int number_of_buckets = T.rows();
	igl::SpatialHash<DerivedV> spatial_hash(s_grid, number_of_buckets);
	for (auto e_it = intersecting_edges.begin(); e_it != intersecting_edges.end(); e_it++) {
		spatial_hash.insert(*e_it);
	}
	for (int f_i = 0; f_i < SF.rows(); f_i++) {
		spatial_hash.insert(Triangle(SF(f_i, 0), SF(f_i, 1), SF(f_i, 2), f_i));
	}

	return spatial_hash.edge_surface_collision(
		intersection_points,
		surface_normals,
		border_points);
}


template< typename DerivedV>
IGL_INLINE bool igl::penetration_depth_and_direction(
	const Eigen::PlainObjectBase< DerivedV > & V,
	const std::vector<Eigen::Vector3d> & intersection_points,
	const std::vector<Eigen::Vector3d> & surface_normals,
	const std::map<int, std::set<int>> & neighbors,
	const std::map<int, std::set<int>> & border_points,
	std::map<int, Eigen::Vector3d> & penetration_direction,
	std::map<int, double> & penetration_depth
)
{
	using namespace Eigen;
	int idx = 0;
	for (auto b_it = border_points.begin(); b_it != border_points.end(); b_it++) {
		const Vector3d p = V.row(b_it->first).transpose();
		double d = 0.;  // penetration depth
		double weight_sum = 0.;
		Vector3d r_tilde = Vector3d::Zero();
		for (auto i_it = b_it->second.begin(); i_it != b_it->second.end(); i_it++) {
			const Vector3d xi = intersection_points[*i_it].transpose();
			const Vector3d ni = surface_normals[*i_it].transpose();
			double weight = 1. / (xi - p).squaredNorm();
			if (isinf(weight)) {
				continue;
			}
			d += weight * (xi - p).dot(ni);
			r_tilde += weight * ni;
			weight_sum += weight;
		}

		if (d != d) {
			std::cout << "weight_sum = " << weight_sum << std::endl;
		}
		d /= weight_sum;
		r_tilde /= weight_sum;
		penetration_depth[b_it->first] = d;
		penetration_direction[b_it->first] = r_tilde.normalized();
		idx++;
	}
	return true;
}

template< typename DerivedV >
IGL_INLINE bool igl::propagation(
	const Eigen::PlainObjectBase< DerivedV > & V,
	const std::map<int, Eigen::Vector3d> & border_penetration_direction,
	const std::map<int, double> & border_penetration_depth,
	const std::set<int> & colliding_points,
	const std::map<int, std::set<int>> & border_points,
	const std::map<int, std::set<int>> & neighbors,
	std::map<int, Eigen::Vector3d> & penetration_direction,
	std::map<int, double> & penetration_depth
)
{
	using namespace Eigen;
	// border_points holds the processed point index
	std::set<int> unprocessed_points = colliding_points;
	std::set<int> border_points_set;
	for (auto it = border_points.begin(); it != border_points.end(); it++) {
		unprocessed_points.erase(it->first);
		border_points_set.insert(it->first);
		penetration_direction[it->first] = border_penetration_direction.at(it->first);
		penetration_depth[it->first] = border_penetration_depth.at(it->first);
	}
	while (!unprocessed_points.empty()) {
		// identify new border points
		std::map<int, std::set<int>> new_border_points; // map key: new border points, second: adjacent processed points
		for (auto it = unprocessed_points.begin(); it != unprocessed_points.end(); it++) {
			for (auto n_it = neighbors.at(*it).begin(); n_it != neighbors.at(*it).end(); n_it++) {
				if (border_points_set.find(*n_it) != border_points_set.end()) { // new border points
					new_border_points[*it].insert(*n_it);
				}
			}
		}
		if (new_border_points.empty()) {
			return true;
		}
		for (auto it = new_border_points.begin(); it != new_border_points.end(); it++) {
			const Vector3d p = V.row(it->first).transpose();
			double d = 0.;
			double weight_sum = 0.;
			Vector3d r_tilde = Vector3d::Zero();
			for (auto p_it = it->second.begin(); p_it != it->second.end(); p_it++) {
				const Vector3d pj = V.row(*p_it).transpose();
				const Vector3d rj = penetration_direction.at(*p_it);
				double dj = penetration_depth.at(*p_it);
				double weight = 1. / (pj - p).squaredNorm();
				if (isinf(weight)) {
					continue;
				}
				d += weight * (pj - p).dot(rj) + dj;
				r_tilde += weight * rj;
				weight_sum += weight;
			}

			d /= weight_sum;
			r_tilde /= weight_sum;

			penetration_direction[it->first] = r_tilde;
			penetration_depth[it->first] = d;
			unprocessed_points.erase(it->first);
			border_points_set.insert(it->first);
		}
	}
	return true;
}




#ifdef IGL_STATIC_LIBRARY
template bool igl::self_collision< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::Matrix<int, -1, 1, 0, -1, 1> const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::vector< Eigen::Triplet<double, int>, std::allocator< Eigen::Triplet<double, int> > > &);
template bool igl::self_collision< Eigen::Matrix<double, -1, -1, 0, -1, -1>,  Eigen::Matrix<int, -1, -1, 0, -1, -1>,  Eigen::Matrix<int, -1, -1, 0, -1, -1> >( Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,  Eigen::Matrix<int, -1, 1, 0, -1, 1> const &,  Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,  Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,  std::set<int, struct std::less<int>,  std::allocator<int> > const &,  Eigen::Matrix<int, -1, 1, 0, -1, 1> &,  Eigen::Matrix<double, -1, 3, 0, -1, 3> &);
template bool igl::point_collisions< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::set<int, struct std::less<int>, std::allocator<int> > const &, std::set<int, struct std::less<int>, std::allocator<int> > &);
template bool igl::self_collision< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::map<int, std::set<int, struct std::less<int>, std::allocator<int> >, struct std::less<int>, std::allocator<struct std::pair<int const, std::set<int, struct std::less<int>, std::allocator<int> > > > > const &, std::set<int, struct std::less<int>, std::allocator<int> > const &, struct igl::RBCData &);



#endif