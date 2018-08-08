#ifndef IGL_STAIRS_H
#define IGL_STAIRS_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
	struct Stairs {
	bool enabled = false;
	int number_of_stairs = 10;
	double step_width = 1.5;
	double step_height = 1.0;
	double start_height = 0.0;
	double start_width = 2.0;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	Stairs() {
		V.resize(4 * number_of_stairs, 3);
		F.resize(2 * number_of_stairs, 3);
		for (int i = 0; i < number_of_stairs; i++) {
			V.block(i * 4, 0, 4, 3) << -10., (double)(-i * step_height - start_height), (double)((i + 1)* step_width - start_width),
				10., (double)(-i * step_height - start_height), (double)((i + 1) * step_width - start_width),
				10., (double)(-i * step_height - start_height), (double)(i * step_width - start_width),
				-10., (double)(-i * step_height - start_height),(double)(i * step_width - start_width);

			F.block(i * 2, 0, 2, 3) << 0 + i * 4, 1 + i * 4, 3 + i * 4,
										1 + i * 4, 2 + i * 4, 3 + i * 4;

		}
	}
};
}

#ifndef IGL_STATIC_LIBRARY
#  include "stairs.cpp"
#endif
#endif

