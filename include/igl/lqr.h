#ifndef IGL_LQR_H
#define IGL_LQR_H
#include "igl_inline.h"
#include <vector>
#include <map>

namespace igl
{	
	// targets key: frame -> val : 3d coordinate
	// frame starts from 0.
	template <
		typename DerivedX,
		typename Scalar>
	IGL_INLINE void lqr(
		int N,
		Scalar dt,
		Scalar rho,
		const DerivedX & x_init,
		const DerivedX & x_target,
		std::vector<DerivedX> &x_sol);
}

#ifndef IGL_STATIC_LIBRARY
#  include "lqr.cpp"
#endif

#endif
