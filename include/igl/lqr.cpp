#include "lqr.h"

#include <vector>
#include <Eigen/Dense>

template <typename DerivedV>
IGL_INLINE void igl::lqr(
	const std::map<int, DerivedV> &x_,
	std::vector<DerivedV> &x)
{
	using namespace std;
	using namespace Eigen;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
#endif
