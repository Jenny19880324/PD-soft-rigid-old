#include "self_collision.h"
#include <Eigen/Geometry>
#include <igl/per_face_normals.h>
#include <iostream>


template<typename DerivedV, typename DerivedF, typename Scalar>
IGL_INLINE bool igl::self_collision(
	const Eigen::PlainObjectBase <DerivedV> & V,
	const Eigen::PlainObjectBase <DerivedF> & F,
	const Eigen::VectorXi & SV,
	const Eigen::PlainObjectBase < DerivedF > & SF,
	std::vector<Eigen::Triplet<Scalar>> & L_triplets)
{
	using namespace std;
	using namespace Eigen;
	typedef Matrix<typename DerivedV::Scalar, 1, 3> RowVectorV3;

	MatrixXd N;
	per_face_normals_stable(V, SF, N);

	for (int i = 0; i < SV.rows(); i++) {
		const RowVectorV3 vtx = V.row(SV(i));
		for (int f_i = 0; f_i < SF.rows(); f_i++) {
			const RowVectorV3 p0 = V.row(SF(f_i, 0));
			const RowVectorV3 r = vtx - p0;
			if (r.dot(N(f_i)) < 0.0) {
				std::cout << "vertex " << SV(i) << "collides with facet " << f_i << std::endl;
			}

		}
	}
}

#ifdef IGL_STATIC_LIBRARY


#endif
