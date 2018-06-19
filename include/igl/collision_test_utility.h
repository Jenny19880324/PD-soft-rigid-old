
#ifndef IGL_COLLISION_TEST_UTILITY_H
#define IGL_COLLISION_TEST_UTILITY_H


#include "igl_inline.h"
#include <igl/rbc.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <array>
#include <unordered_set>


namespace igl
{
	typedef int Point;

	struct Edge {

		Edge(int _from, int _to, int _idx) : from(_from), to(_to), idx(_idx) {}
		Edge(const Edge & other) { from = other.from; to = other.to; idx = other.idx; }
		Edge &operator=(const Edge &other) { from = other.from; to = other.to; idx = other.idx; return *this; }
		bool operator==(const Edge &other) const {
			return idx == other.idx;
		}
		// compare for order
		bool operator < (const Edge & e) const
		{
			return (from < e.from) || ((!(e.from < from)) && (to < e.to));
		}

		int from, to;
		int idx;
	};

	struct Triangle {
		Triangle() {}
		Triangle(int _a, int _b, int _c, int _idx) : a(_a), b(_b), c(_c), idx(_idx) {}
		Triangle(const Triangle &other) { a = other.a; b = other.b; c = other.c; idx = other.idx; }
		Triangle &operator=(const Triangle &other) { a = other.a; b = other.b; c = other.c; idx = other.idx; return *this; }
		bool operator==(const Triangle &other) const {
			return idx == other.idx;
		}
		int a, b, c;
		int idx;
	};

	struct Tetrahedron {
		Tetrahedron(int _a, int _b, int _c, int _d, int _idx) : a(_a), b(_b), c(_c), d(_d), idx(_idx) {}
		Tetrahedron(const Tetrahedron & other) { a = other.a; b = other.b; c = other.c; d = other.d; idx = other.idx; }
		Tetrahedron &operator=(const Tetrahedron &other) { a = other.a; b = other.b; c = other.c; d = other.d; idx = other.idx; return *this; }
		bool operator==(const Tetrahedron &other) const {
			return idx == other.idx;
		}
		int a, b, c, d;
		int idx;
	};

	struct CellCoordinate {
		CellCoordinate(int _i, int _j, int _k) : i(_i), j(_j), k(_k) {}
		bool operator==(const CellCoordinate &other) const {
			return (i == other.i &&
				j == other.j &&
				k == other.k);
		}
		int i, j, k;
	};


	struct EdgeHasher
	{
		size_t operator()(const Edge & e) const
		{
			return e.idx;
		}
	};

	struct TriangleHasher
	{
		size_t operator()(const Triangle &tri) const
		{
			return tri.idx;
		}
	};

	struct TetrahedronHasher
	{
		size_t operator()(const Tetrahedron &tet) const
		{
			return tet.idx;
		}
	};


	struct Primitives {
		std::unordered_set<Point>                          point_set;
		std::unordered_set<Edge, EdgeHasher>               edge_set;
		std::unordered_set<Triangle, TriangleHasher>       triangle_set;
		std::unordered_set<Tetrahedron, TetrahedronHasher> tetrahedron_set;
	};





	template<typename DerivedT>
	IGL_INLINE void find_neighbors(
		const Eigen::PlainObjectBase< DerivedT > & T,
		std::map<int, std::set<int>> & neighbors);

	template<typename Scalar>
	IGL_INLINE bool point_in_tet(
		const Eigen::RowVector3d & vtx,
		const Eigen::Matrix<Scalar, 4, 3> & tet);

	template<typename Scalar>
	IGL_INLINE bool edge_surface_intersection(
		const Eigen::Matrix<Scalar, 2, 3> & edge,
		const Eigen::Matrix<Scalar, 3, 3> & tri,
		Scalar &t);

	template<typename DerivedV, typename DerivedT>
	IGL_INLINE double compute_grid_size(
		const Eigen::PlainObjectBase< DerivedV > & V,
		const Eigen::PlainObjectBase< DerivedT > & T
	);

}
#ifndef IGL_STATIC_LIBRARY
#   include "self_collision.cpp"
#endif
#endif //IGL_SEGMENT_SEGMENT_INTERSECT_H
