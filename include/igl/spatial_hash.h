#ifndef IGL_SPATIAL_HASH_H
#define IGL_SPATIAL_HASH_H

#include <list>
#include <algorithm>
#include <functional>
#include <array>
#include <Eigen/Core>
#include <math.h>
#include <igl/collision_test_utility.h>

namespace igl
{
template <typename DerivedV>
	class SpatialHash {
		private:
			struct SpatialHashEntry {
				SpatialHashEntry(const CellCoordinate &_key, const Primitives & _val):key(_key), val(_val){}
				CellCoordinate key;
				Primitives val;
			};
			
			struct AABB {
				double minX, minY, minZ;
				double maxX, maxY, maxZ;
				
				AABB(const Triangle       &tri, const Eigen::MatrixBase<DerivedV> & V);
				AABB(const Tetrahedron    &tet, const Eigen::MatrixBase<DerivedV> & V);
				AABB(const CellCoordinate &cc,  const double s_grid);
				
				IGL_INLINE bool intersect(const AABB &other);
			};
			
			int H(int i, int j, int k) { return (i * 3117209) ^ (j * 6656291) ^ (k * 15485783);}

			std::vector<std::list<SpatialHashEntry>> table;
			double s_grid;
			int number_of_buckets;
		public:

	
			SpatialHash(const double _s_grid, const int _number_of_buckets):
				s_grid(_s_grid), 
				number_of_buckets(_number_of_buckets) {
				table.resize(number_of_buckets);
			}
			
			IGL_INLINE bool insert(const Point       & p);
			IGL_INLINE bool insert(const Edge        & e);
			IGL_INLINE bool insert(const Triangle    & tri);
			IGL_INLINE bool insert(const Tetrahedron & tet);

			IGL_INLINE bool point_tet_collision(
				const std::set<int> & b_set,
				std::set<int> & colliding_points);
			IGL_INLINE bool edge_surface_collision(
				std::vector<Eigen::Vector3d> & intersection_points,
				std::vector<Eigen::Vector3d> & surface_normals,
				std::map<int, std::set<int>> & border_points);
	};
	
	
}


#ifndef IGL_STATIC_LIBRARY
#  include "spatial_hash.cpp"
#endif

#endif
