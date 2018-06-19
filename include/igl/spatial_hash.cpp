#include "spatial_hash.h"
#include "collision_test_utility.h"
#include <cmath>

extern Eigen::MatrixXd U;

template < typename DerivedV >
igl::SpatialHash<DerivedV>::AABB::AABB(const Triangle &tri, const Eigen::MatrixBase<DerivedV> & V)
{
	Eigen::Matrix3d tri_V;
	tri_V << U.row(tri.a),
			 U.row(tri.b),
			 U.row(tri.c);
	Eigen::Vector3d m = tri_V.colwise().minCoeff();
	Eigen::Vector3d M = tri_V.colwise().maxCoeff();
	
	minX = m.x();
	minY = m.y();
	minZ = m.z();
	maxX = M.x();
	maxY = M.y();
	maxZ = M.z();
}

template < typename DerivedV >
igl::SpatialHash<DerivedV>::AABB::AABB(const Tetrahedron &tet, const Eigen::MatrixBase<DerivedV> & V)
{
	Eigen::Matrix<double, 4, 3> tet_V;
	tet_V << U.row(tet.a),
			 U.row(tet.b),
			 U.row(tet.c),
			 U.row(tet.d);
	Eigen::Vector3d m = tet_V.colwise().minCoeff();
	Eigen::Vector3d M = tet_V.colwise().maxCoeff();
	
	minX = m.x();
	minY = m.y();
	minZ = m.z();
	maxX = M.x();
	maxY = M.y();
	maxZ = M.z();
}

template< typename DerivedV >
igl::SpatialHash<DerivedV>::AABB::AABB(const CellCoordinate &cc, const double s_grid)
{
	minX = cc.i * s_grid;
	minY = cc.j * s_grid;
	minZ = cc.k * s_grid;
	maxX = (cc.i + 1) * s_grid;
	maxY = (cc.j + 1) * s_grid;
	maxZ = (cc.k + 1) * s_grid;
}

template< typename DerivedV >
IGL_INLINE bool igl::SpatialHash<DerivedV>::AABB::intersect(const AABB &other)
{
	return  (minX < other.maxX && maxX > other.minX) &&
			(minY < other.maxY && maxY > other.minY) &&
			(minZ < other.maxZ && maxZ > other.minZ);
}


template < typename DerivedV >
IGL_INLINE bool igl::SpatialHash<DerivedV>::insert(const Point &p) 
{
	int x = (int)floor(U.row(p).x() / s_grid);
	int y = (int)floor(U.row(p).y() / s_grid);
	int z = (int)floor(U.row(p).z() / s_grid);
	int h = (H(x, y, z) % number_of_buckets + number_of_buckets) % number_of_buckets;
	CellCoordinate key = CellCoordinate(x, y, z);
	
	
	if (table[h].empty()) {
		Primitives val;
		val.point_set.insert(p);
		table[h].push_back(SpatialHashEntry(key, val));
		return true;
	}
	else {
		for (auto it = table[h].begin(); it != table[h].end(); it++)
		{
			if (it->key == key) {
				auto res = it->val.point_set.insert(p);
				if (res.second) {
					return true;
				}
				else {
					return false;
				}
			}
		}
		Primitives val;
		val.point_set.insert(p);
		table[h].push_back(SpatialHashEntry(key, val));
	}
}


template< typename DerivedV >
IGL_INLINE bool igl::SpatialHash<DerivedV>::insert(const Edge &e)
{
	int xs = (int)floor(U.row(e.from).x() / s_grid);
	int ys = (int)floor(U.row(e.from).y() / s_grid);
	int zs = (int)floor(U.row(e.from).z() / s_grid);
	int xe = (int)floor(U.row(e.to).x()   / s_grid);
	int ye = (int)floor(U.row(e.to).y()   / s_grid);
	int ze = (int)floor(U.row(e.to).z()   / s_grid);
	
	int x = floor(xs);
	int y = floor(ys);
	int z = floor(zs);
	
	double dx = xe - xs;
	double dy = ye - ys;
	double dz = ze - zs;
	
	double tx = (xs - x) * dx;
	double ty = (ys - y) * dy;
	double tz = (zs - z) * dz;
	
	int n = abs(floor(xe) - x) + abs(floor(ye) - y) + abs(floor(ze) - z);
	
	for (int i = 0; i < n; i++) {
		// insert edge under key (x, y, z)
		assert(x <= xe);
		assert(y <= ye);
		assert(z <= ze);
		int h = (H(x, y, z) % number_of_buckets + number_of_buckets) % number_of_buckets;
		CellCoordinate key = CellCoordinate(x, y, z);
		assert(h >= 0 && h < number_of_buckets);
		if (table[h].empty()) {
			Primitives val;
			val.edge_set.insert(e);
			table[h].push_back(SpatialHashEntry(key, val));
		}
		else {
			bool key_found = false;
			for (auto it = table[h].begin(); it != table[h].end(); it++) {
				if (it->key == key) {
					key_found = true;
					auto res = it->val.edge_set.insert(e);
					if (!res.second) {
						//std::cout << "edge already inserted!" << std::endl;
						return false;
					}
				}
			}
			if (!key_found) {
				//std::cout << "insert edge hash collision" << std::endl;
				Primitives val;
				val.edge_set.insert(e);
				table[h].push_back(SpatialHashEntry(key, val));
			}
		}
		
		// incremental traversal
		if (tx < ty) {
			if (tx < tz) {
				x = x + ( dx > 0 ? 1 : -1);
				tx = tx + abs(1 / dx);
			}
			else {
				z = z + (dz > 0 ? 1 : -1);
				tz = tz + abs(1 / dz);
			}
		}
		else {
			if (ty < tz) {
				y = y + (dy > 0 ? 1 : -1);
				ty = ty + abs(1 / dy);
			}
			else {
				z = z + (dz > 0 ? 1 : -1);
				tz = tz + abs(1 / dz);
			}
		}
	}
	return true;
}


template < typename DerivedV >
IGL_INLINE bool igl::SpatialHash<DerivedV>::insert(const Triangle & tri)
{
	AABB box = AABB(tri, U);
	int xs = (int)floor(box.minX / s_grid);
	int ys = (int)floor(box.minY / s_grid);
	int zs = (int)floor(box.minZ / s_grid);
	int xe = (int)floor(box.maxX / s_grid);
	int ye = (int)floor(box.maxY / s_grid);
	int ze = (int)floor(box.maxZ / s_grid);
	
	Eigen::RowVector3d e1 = U.row(tri.a) - U.row(tri.c);
	Eigen::RowVector3d e2 = U.row(tri.b) - U.row(tri.c);
	Eigen::RowVector3d n = e1.cross(e2);
	n.normalize();
	double xn, yn, zn;
	double xp, yp, zp;
	double x0 = U.row(tri.a).x();
	double y0 = U.row(tri.b).y();
	double z0 = U.row(tri.c).z();
	double A = n.x();
	double B = n.y();
	double C = n.z();
	double D = - n.x() * x0 - n.y() * y0 - n.z() * z0;
	
	
	for (int x = xs; x <= xe; x++) {
		for (int y = ys; y <= ye; y++) {
			for (int z = zs; z <= ze; z++) { // for each cell that is intersected by the AABB
				// test if the same cell is intersected by the plane defined by the triangle
				// box-plane intersection test [Gre94]
				if (n.x() > 0.) {
					xp = (x + 1) * s_grid;
					xn = x * s_grid;
				}
				else {
					xp = x * s_grid;
					xn = (x + 1) * s_grid;
				}

				if (n.y() > 0.) {
					yp = (y + 1) * s_grid;
					yn = y * s_grid;
				}
				else {
					yp = y * s_grid;
					yn = (y + 1) * s_grid;
				}

				if (n.z() > 0.) {
					zp = (z + 1) * s_grid;
					zn = z * s_grid;
				}
				else {
					zp = z * s_grid;
					zn = (z + 1) * s_grid;
				}
				
				if((A * xp + B * yp + C * zp + D > 0) && (A * xn + B * yn + C * zn + D < 0)){
					int h = (H(x, y, z) % number_of_buckets + number_of_buckets) % number_of_buckets;
					CellCoordinate key = CellCoordinate(x, y, z);
					assert(h >= 0 && h < number_of_buckets);
					if (table[h].empty()) {
						Primitives val;
						val.triangle_set.insert(tri);
						table[h].push_back(SpatialHashEntry(key, val));
					}
					else{
						bool key_found = false;
						for (auto it = table[h].begin(); it != table[h].end(); it++) {
							if (it->key == key) {
								key_found = true;
								auto res = it->val.triangle_set.insert(tri);
								if (!res.second) {
									//std::cout << "triangle already inserted!" << std::endl;
									return false;
								}
							}
						}
						if (!key_found) {
							//std::cout << "insert triangle hash collision" << std::endl;
							Primitives val;
							val.triangle_set.insert(tri);
							table[h].push_back(SpatialHashEntry(key, val));
						}
					}
				}
			}
		}
	}
	return true;

}


template < typename DerivedV >
IGL_INLINE bool igl::SpatialHash<DerivedV>::insert(const Tetrahedron & tet)
{
	AABB box = AABB(tet, U);
	int xs = (int)floor(box.minX / s_grid);
	int ys = (int)floor(box.minY / s_grid);
	int zs = (int)floor(box.minZ / s_grid);
	int xe = (int)floor(box.maxX / s_grid);
	int ye = (int)floor(box.maxY / s_grid);
	int ze = (int)floor(box.maxZ / s_grid);
	
	for (int x = xs; x <= xe; x++) {
		for (int y = ys; y <= ye; y++) {
			for (int z = zs; z <= ze; z++) {
				int h = (H(x, y, z) % number_of_buckets + number_of_buckets) % number_of_buckets;
				CellCoordinate key = CellCoordinate(x, y, z);
				assert(h >= 0 && h < number_of_buckets);
				if (table[h].empty()) {
					Primitives val;
					val.tetrahedron_set.insert(tet);
					table[h].push_back(SpatialHashEntry(key, val));
				}
				else {
					bool key_found = false;
					for (auto it = table[h].begin(); it != table[h].end(); it++) {
						if (it->key == key) {
							key_found = true;
							auto res = it->val.tetrahedron_set.insert(tet);
							if (!res.second) {
								//std::cout << "tetrahedron already inserted!" << std::endl;
								return false;
							}
						}
					}
					if (!key_found) {
						Primitives val;
						val.tetrahedron_set.insert(tet);
						table[h].push_back(SpatialHashEntry(key, val));
					}
				}
			}
		}
	}
	
	return true;
	
}


template < typename DerivedV >
IGL_INLINE bool igl::SpatialHash<DerivedV>::point_tet_collision(
	const std::set<int> & b_set,
	std::set<int> & colliding_points
)
{
	bool is_intersected = false;
	for (int tb_i = 0; tb_i < table.size(); tb_i++) {
		for (auto it = table[tb_i].begin(); it != table[tb_i].end(); it++) {
			for (auto p_it = it->val.point_set.begin(); p_it != it->val.point_set.end(); p_it++) {
				for (auto tet_it = it->val.tetrahedron_set.begin(); tet_it != it->val.tetrahedron_set.end(); tet_it++) {
					if (b_set.find(*p_it) != b_set.end()) {
						continue;
					}
					const Eigen::RowVector3d vtx = U.row(*p_it);
					Eigen::Matrix<double, 4, 3> tet;
					tet << U.row((*tet_it).a),
						U.row((*tet_it).b),
						U.row((*tet_it).c),
						U.row((*tet_it).d);

					if (point_in_tet(vtx, tet)) {
						colliding_points.insert(*p_it);
						is_intersected = true;
					}
				}
			}
		}
	}

	return is_intersected;
}


template < typename DerivedV >
IGL_INLINE bool igl::SpatialHash<DerivedV>::edge_surface_collision(
	std::vector<Eigen::Vector3d> & intersection_points,
	std::vector<Eigen::Vector3d> & surface_normals,
	std::map<int, std::set<int>> & border_points
)
{
	using namespace Eigen;
	bool is_intersected = false;
	for (int tb_i = 0; tb_i < table.size(); tb_i++) {
		for (auto it = table[tb_i].begin(); it != table[tb_i].end(); it++) {
			for (auto e_it = it->val.edge_set.begin(); e_it != it->val.edge_set.end(); e_it++) {
				double max_t = -1;
				Triangle max_tri;
				Eigen::Matrix<double, 2, 3> edge;
				edge << U.row(e_it->from),
						U.row(e_it->to);
				for (auto s_it = it->val.triangle_set.begin(); s_it != it->val.triangle_set.end(); s_it++) {


					Eigen::Matrix<double, 3, 3> tri;
					double t = -1;
					tri << U.row((*s_it).a),
						U.row((*s_it).b),
						U.row((*s_it).c);

					if (edge_surface_intersection(edge, tri, t)) {
						if (t > max_t) {
							max_t = t;
							max_tri = *s_it;
						}
						is_intersected = true;
					}
				}
				
				if (max_t <= 0 || max_t >= 1) {
					continue;
				}

				// assuming all normals of SF is pointed outward
				RowVector3d p = U.row(e_it->from);
				RowVector3d p_end = U.row(e_it->to);
				RowVector3d d = p_end - p;
				RowVector3d intersection_point = p + max_t * d;
				//std::cout << "max_t = " << max_t << std::endl;

				Vector3d r1 = (U.row(max_tri.b) - U.row(max_tri.a)).transpose();
				Vector3d r2 = (U.row(max_tri.c) - U.row(max_tri.a)).transpose();
				Vector3d surface_normal = (r1.cross(r2)).normalized();

				border_points[e_it->from].insert(intersection_points.size());
				intersection_points.push_back(intersection_point);
				surface_normals.push_back(surface_normal);
			}
		}
	}
	
	return is_intersected;
}
#ifdef IGL_STATIC_LIBRARY
template bool igl::SpatialHash< Eigen::Matrix<double, -1, -1, 0, -1, -1> >::insert(int const &);
template bool igl::SpatialHash< Eigen::Matrix<double, -1, -1, 0, -1, -1> >::insert(struct igl::Edge const &);
template bool igl::SpatialHash< Eigen::Matrix<double, -1, -1, 0, -1, -1> >::insert(struct igl::Triangle const &);
template bool igl::SpatialHash< Eigen::Matrix<double, -1, -1, 0, -1, -1> >::insert(struct igl::Tetrahedron const &);
template bool igl::SpatialHash<Eigen::Matrix<double, -1, -1, 0, -1, -1> >::point_tet_collision(std::set<int, struct std::less<int>, std::allocator<int> > const &, std::set<int, struct std::less<int>, std::allocator<int> > &);
template bool igl::SpatialHash<Eigen::Matrix<double, -1, -1, 0, -1, -1> >::edge_surface_collision(std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, std::allocator<  Eigen::Matrix<double, 3, 1, 0, 3, 1> > > &, std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, std::allocator<Eigen::Matrix<double, 3, 1, 0, 3, 1> > > &, std::map<int, std::set<int, struct std::less<int>, std::allocator<int> >, struct std::less<int>, std::allocator<struct std::pair<int const, std::set<int, struct std::less<int>, std::allocator<int> > > > > &);


#endif