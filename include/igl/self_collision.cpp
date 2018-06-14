#include "self_collision.h"
#include <Eigen/Geometry>
#include <igl/per_face_normals.h>
#include <iostream>
#include <algorithm>


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


template<typename DerivedV, typename DerivedT>
IGL_INLINE bool igl::self_collision(
	const Eigen::PlainObjectBase < DerivedV > & V,
	const Eigen::VectorXi & SV,
	const Eigen::PlainObjectBase < DerivedT > & ST,
	const std::set<int> &b_set,
	Eigen::VectorXi & b,
	Eigen::MatrixX3d & bc)
{
	using namespace std;
	using namespace Eigen;
	typedef Matrix<typename DerivedV::Scalar, 1, 3> RowVectorV3;

	b.resize(0);
	bc.resize(0, Eigen::NoChange);

	for (int i = 0; i < SV.rows(); i++) {
		if (b_set.find(SV(i)) != b_set.end()) {
			continue;
		}
		double x = V(SV(i), 0);
		double y = V(SV(i), 1);
		double z = V(SV(i), 2);
		for (int tet_i = 0; tet_i < ST.rows(); tet_i++) {
			int t1 = ST(tet_i, 0);
			int t2 = ST(tet_i, 1);
			int t3 = ST(tet_i, 2);
			int t4 = ST(tet_i, 3);

			double x1 = V(t1, 0);
			double y1 = V(t1, 1);
			double z1 = V(t1, 2);

			double x2 = V(t2, 0);
			double y2 = V(t2, 1);
			double z2 = V(t2, 2);

			double x3 = V(t3, 0);
			double y3 = V(t3, 1);
			double z3 = V(t3, 2);

			double x4 = V(t4, 0);
			double y4 = V(t4, 1);
			double z4 = V(t4, 2);


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

			if (D0.determinant() > 0.) {
				std::cout << "tet_i = " << tet_i << std::endl;
			}
			double eps = 1e-8;
			double det[5] = { D0.determinant(), D1.determinant(), D2.determinant(), D3.determinant(), D4.determinant() };
			if ((det[0] > eps && det[1] > eps && det[2] > eps && det[3] > eps && det[4] > eps) ||
				(det[0] < -eps && det[1] < -eps && det[2] < -eps && det[3] < -eps && det[4] < -eps)) {
				std::cout << "vertex " << SV(i) << " collides with tet of vert " << t1 << ", " << t2 << ", " << t3 << ", " << t4 << std::endl;
				assert(det[0] = det[1] + det[2] + det[3] + det[4]);
				double min_h = 1e8;
				Vector3d min_p1, min_p2, min_p3;
				for (int v_i = 0; v_i < 4; v_i++) {
					const Vector3d p1 = D0.block((v_i + 1) % 4, 0, 1, 3).transpose();
					const Vector3d p2 = D0.block((v_i + 2) % 4, 0, 1, 3).transpose();
					const Vector3d p3 = D0.block((v_i + 3) % 4, 0, 1, 3).transpose();

					double A = fabs((0.5 * (p3 - p1).cross(p3 - p2)).norm());
					double h = 3 * fabs(det[v_i + 1]) / A;

					if (h < min_h) {
						min_h = h;
						min_p1 = p1;
						min_p2 = p2;
						min_p3 = p3;
					}
				}

				// project a point to the triangle
				Vector3d p = Vector3d(x, y, z);
				Vector3d u = min_p2 - min_p1;
				Vector3d v = min_p3 - min_p1;
				Vector3d n = u.cross(v);
				Vector3d w = p - min_p1;
				Vector3d p_projected = p;

				double gamma = (u.cross(w)).dot(n) / (n.squaredNorm());
				double beta = (w.cross(v)).dot(n) / (n.squaredNorm());
				double alpha = 1. - gamma - beta;

				if (alpha >= 0 && alpha <= 1 &&
					beta >= 0 && beta <= 1 &&
					gamma >= 0 && gamma <= 1) {
					std::cout << "The point p lies inside T" << std::endl;
					p_projected = alpha * min_p1 + beta * min_p2 + gamma * min_p3;
					b.conservativeResize(b.rows() + 1);
					b(b.rows() - 1) = SV(i);
					bc.conservativeResize(b.rows(), Eigen::NoChange);
					bc.row(b.rows() - 1) = p_projected.transpose();
					//bc.row(b.rows() - 1) = p.transpose();

					std::cout << "min_h = " << min_h << std::endl;
					std::cout << "b = " << i << std::endl;
					std::cout << "bc = " << p_projected << std::endl;
					std::cout << "p = " << p << std::endl;
				}
			}
		}
	}
	return false;
}

#ifdef IGL_STATIC_LIBRARY
template bool igl::self_collision< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::Matrix<int, -1, 1, 0, -1, 1> const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::vector< Eigen::Triplet<double, int>, std::allocator< Eigen::Triplet<double, int> > > &);
template bool igl::self_collision< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::Matrix<int, -1, 1, 0, -1, 1> const &, Eigen::PlainObjectBase< Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, std::set<int, struct std::less<int>, std::allocator<int> > const &, Eigen::Matrix<int, -1, 1, 0, -1, 1> &, Eigen::Matrix<double, -1, 3, 0, -1, 3> &);
#endif