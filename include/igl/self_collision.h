
#ifndef IGL_SELF_COLLISION_H
#define IGL_SELF_COLLISION_H


#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
namespace igl
{
	template<typename DerivedV, typename DerivedT, typename Scalar>
	IGL_INLINE bool self_collision(
		const Eigen::PlainObjectBase <DerivedV> & V,
		const Eigen::VectorXi & SV,
		const Eigen::PlainObjectBase <DerivedT> & ST,
		std::vector<Eigen::Triplet<Scalar>> & L_triplets);


	template<typename DerivedV, typename DerivedT>
	IGL_INLINE bool self_collision(
		const Eigen::PlainObjectBase < DerivedV > & V,
		const Eigen::VectorXi & SV,
		const Eigen::PlainObjectBase < DerivedT > & ST,
		const std::set<int> &b_set,
		Eigen::VectorXi & b,
		Eigen::MatrixX3d & bc);


}
#ifndef IGL_STATIC_LIBRARY
#   include "self_collision.cpp"
#endif
#endif //IGL_SEGMENT_SEGMENT_INTERSECT_H
