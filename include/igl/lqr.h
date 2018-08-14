#ifndef IGL_LQR_H
#define IGL_LQR_H
#include "igl_inline.h"
#include <vector>
#include <map>

namespace igl
{	
	template <typename DerivedV>
	IGL_INLINE void lqr(
		const std::map<int, DerivedV> &x_,
		std::vector<DerivedV> &x);
}

#ifndef IGL_STATIC_LIBRARY
#  include "lqr.cpp"
#endif

#endif
