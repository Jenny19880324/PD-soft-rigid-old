
#ifndef IGL_SELF_COLLISION_H
#define IGL_SELF_COLLISION_H


#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
namespace igl
{
	template<typename DerivedV, typename DerivedF, typename Scalar>
	IGL_INLINE bool self_collision(
		const Eigen::PlainObjectBase <DerivedV> & V,
		const Eigen::PlainObjectBase <DerivedF> & F,
		const Eigen::VectorXi & SV,
		const Eigen::PlainObjectBase <DerivedF> & SF,
		std::vector<Eigen::Triplet<Scalar>> & L_triplets);

}
#ifndef IGL_STATIC_LIBRARY
#   include "self_collision.cpp"
#endif
#endif //IGL_SEGMENT_SEGMENT_INTERSECT_H
