
#ifndef IGL_SELF_COLLISION_H
#define IGL_SELF_COLLISION_H


#include "igl_inline.h"
#include <igl/rbc.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include "collision_test_utility.h"
#include "spatial_hash.h"


namespace igl
{

	template<typename DerivedV, typename DerivedT, typename Scalar>
	IGL_INLINE bool self_collision(
		const Eigen::PlainObjectBase <DerivedV> & V,
		const Eigen::VectorXi & SV,
		const Eigen::PlainObjectBase <DerivedT> & ST,
		std::vector<Eigen::Triplet<Scalar>> & L_triplets);


	template<typename DerivedV, typename DerivedT, typename DerivedF>
	IGL_INLINE bool self_collision(
		const Eigen::PlainObjectBase < DerivedV > & V,
		const Eigen::VectorXi & SV,
		const Eigen::PlainObjectBase < DerivedT > & ST,
		const Eigen::PlainObjectBase < DerivedF > & SF,
		const std::set<int> &b_set,
		Eigen::VectorXi & b,
		Eigen::MatrixX3d & bc);


	template<typename DerivedV, typename DerivedT, typename DerivedF>
	IGL_INLINE bool self_collision(
		const Eigen::PlainObjectBase< DerivedV > & V,
		const Eigen::PlainObjectBase< DerivedT > & T,
		const Eigen::PlainObjectBase< DerivedF > & SF,
		const std::map<int, std::set<int>> & neighbors,
		const std::set<int> & b_set,
		RBCData & data);


	// step 1
	template<typename DerivedV, typename DerivedT>
	IGL_INLINE bool point_collisions(
		const Eigen::PlainObjectBase< DerivedV > & V,
		const Eigen::PlainObjectBase< DerivedT > & T,
		const std::set<int> & b_set,
		std::set<int> & colliding_points
	);

	// step 2
	template<typename DerivedV, typename DerivedT, typename DerivedF>
	IGL_INLINE bool edge_intersections(
		const Eigen::PlainObjectBase< DerivedV > & V,
		const Eigen::PlainObjectBase< DerivedT > & T,
		const Eigen::PlainObjectBase< DerivedF > & SF,
		const std::map<int, std::set<int>> & neighbors,
		const std::set<int> & colliding_points,
		const std::set<int> & b_set,
		std::vector<Eigen::Vector3d> & intersection_points,
		std::vector<Eigen::Vector3d> & surface_normals,
		std::map<int, std::set<int>> & border_points); //map key: border_point index, second: index in intersection_points and surface_normals

	// step 3
	template<typename DerivedV>
	IGL_INLINE bool penetration_depth_and_direction(
		const Eigen::PlainObjectBase< DerivedV > & V,
		const std::vector<Eigen::Vector3d> & intersection_points,
		const std::vector<Eigen::Vector3d> & surface_normals,
		const std::map<int, std::set<int>> & neighbors,
		const std::map<int, std::set<int>> & border_points,
		std::map<int, Eigen::Vector3d> & penetration_direction,
		std::map<int, double> & penetration_depth
	);

	// step 4
	template< typename DerivedV> 
	IGL_INLINE bool propagation(
		const Eigen::PlainObjectBase< DerivedV > & V,
		const std::map<int, Eigen::Vector3d> & border_penetration_direction,
		const std::map<int, double> & border_penetration_depth,
		const std::set<int> & colliding_points,
		const std::map<int, std::set<int>> & border_points,
		const std::map<int, std::set<int>> & neighbors,
		std::map<int, Eigen::Vector3d> & penetration_direction,
		std::map<int, double> & penetration_depth
	);

}
#ifndef IGL_STATIC_LIBRARY
#   include "self_collision.cpp"
#endif
#endif //IGL_SEGMENT_SEGMENT_INTERSECT_H
