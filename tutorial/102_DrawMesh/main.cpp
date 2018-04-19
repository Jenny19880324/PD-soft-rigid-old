#include <igl/readMESH.h>
#include <igl/trackball.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/MeshGL.h>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi T;
Eigen::MatrixXd tV;
Eigen::MatrixXi tF;
Eigen::MatrixXd tC;


struct SlicePlane {
	bool enabled = false;
	bool active = false;
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	Eigen::Vector3d normal = Eigen::Vector3d(0., 0., 1.);
	Eigen::Vector3d axis1  = Eigen::Vector3d(1., 0., 0.);
	Eigen::Vector3d axis2  = Eigen::Vector3d(0., 1., 0.);
	Eigen::Matrix4d trans  = Eigen::Matrix4d::Identity();
	Eigen::MatrixXd vertices;
	Eigen::MatrixXd stored_vertices;

	SlicePlane() {
		vertices.resize(4, 3);
		vertices <<
			-2., 2., 0.,
			2., 2., 0.,
			2., -2., 0.,
			-2., -2., 0.;

		stored_vertices.resize(4, 3);
		stored_vertices <<
			-2., 2., 0.,
			2., 2., 0.,
			2., -2., 0.,
			-2., -2., 0.;
	}

};

SlicePlane slice_plane;

Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
Eigen::Quaternionf down_rotation = Eigen::Quaternionf::Identity();


template <typename DerivedV, typename DerivedT>
void updateTetVisibility(
	const SlicePlane &slice_plane,
	const Eigen::PlainObjectBase<DerivedV> &V,
	const Eigen::PlainObjectBase<DerivedT> &T,
	std::vector<bool> &tet_visibility)
{
	std::vector<bool> vert_visibility(V.rows());
	for (int v_i = 0; v_i < V.rows(); ++v_i) {
		if (slice_plane.enabled) {
			vert_visibility[v_i] = slice_plane.normal.dot(Eigen::Vector3d(V.row(v_i)) - slice_plane.center) > 0.0;
		}
		else {
			vert_visibility[v_i] = true;
		}
	}
	tet_visibility.resize(T.rows());
	for (int tet_i = 0; tet_i < T.rows(); ++tet_i) {
		tet_visibility[tet_i] =
			vert_visibility[T(tet_i, 0)] &&
			vert_visibility[T(tet_i, 1)] &&
			vert_visibility[T(tet_i, 2)] &&
			vert_visibility[T(tet_i, 3)];
	}
}


template <typename DerivedV, typename DerivedT, typename DerivedF, typename DerivedC>
void updateVisibleIndices(
	const SlicePlane &slice_plane,
	const Eigen::PlainObjectBase<DerivedV> &V,
	const Eigen::PlainObjectBase<DerivedT> &T,
	Eigen::PlainObjectBase<DerivedF> &tF,
	Eigen::PlainObjectBase<DerivedC> &tC)
{
	tF.resize(T.rows() * 4, 3); // max out the size for now
	tC.resize(tF.rows(), 4);
	Eigen::MatrixX4d colorScheme;
	colorScheme.resize(12, 4);

	colorScheme <<
		166, 206, 227, 255
		, 31, 120, 180, 255
		, 178, 223, 138, 255
		, 51, 160, 44, 255
		, 251, 154, 153, 255
		, 227, 26, 28, 255
		, 253, 191, 111, 255
		, 255, 127, 0, 255
		, 202, 178, 214, 255
		, 106, 61, 154, 255
		, 255, 255, 153, 255
		, 177, 89, 40, 255;

	colorScheme /= 255.0f;

	std::vector<bool> tet_visibility;
	updateTetVisibility(slice_plane, V, T, tet_visibility);

	size_t tet_visible_idx = 0;

	for (size_t tet_i = 0; tet_i < T.rows(); ++tet_i) {
		if (!tet_visibility[tet_i]) {
			continue;
		}

		tF(4 * tet_visible_idx + 0, 0) = 12 * tet_i + 0;
		tF(4 * tet_visible_idx + 0, 1) = 12 * tet_i + 1;
		tF(4 * tet_visible_idx + 0, 2) = 12 * tet_i + 2;

		tF(4 * tet_visible_idx + 1, 0) = 12 * tet_i + 3;
		tF(4 * tet_visible_idx + 1, 1) = 12 * tet_i + 4;
		tF(4 * tet_visible_idx + 1, 2) = 12 * tet_i + 5;

		tF(4 * tet_visible_idx + 2, 0) = 12 * tet_i + 6;
		tF(4 * tet_visible_idx + 2, 1) = 12 * tet_i + 7;
		tF(4 * tet_visible_idx + 2, 2) = 12 * tet_i + 8;

		tF(4 * tet_visible_idx + 3, 0) = 12 * tet_i + 9;
		tF(4 * tet_visible_idx + 3, 1) = 12 * tet_i + 10;
		tF(4 * tet_visible_idx + 3, 2) = 12 * tet_i + 11;

		const size_t colorIdx = T(tet_i, 4) % 12;
		Eigen::Vector4d color = colorScheme.row(colorIdx);
		tC.row(4 * tet_visible_idx + 0) = color;
		tC.row(4 * tet_visible_idx + 1) = color;
		tC.row(4 * tet_visible_idx + 2) = color;
		tC.row(4 * tet_visible_idx + 3) = color;

		++tet_visible_idx;
	}
	tF.conservativeResize(tet_visible_idx * 4, Eigen::NoChange);
	tC.conservativeResize(tet_visible_idx * 4, Eigen::NoChange);
}


bool mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y) {
	if (slice_plane.active && viewer.down) {
		igl::trackball(
			viewer.core.viewport(2),
			viewer.core.viewport(3),
			2.0f,
			down_rotation,
			viewer.down_mouse_x,
			viewer.down_mouse_y,
			mouse_x,
			mouse_y,
			trackball_angle);

		Eigen::Matrix4f rot = Eigen::Matrix4f::Identity();
		rot.block<3, 3>(0, 0) = trackball_angle.toRotationMatrix();
		slice_plane.trans = (viewer.core.model.inverse() * rot * viewer.core.model).cast<double>();
		slice_plane.axis1 = (slice_plane.trans * Eigen::Vector4d(slice_plane.axis1.x(), slice_plane.axis1.y(), slice_plane.axis1.z(), 1.0)).block<3, 1>(0, 0);
		slice_plane.axis2 = (slice_plane.trans * Eigen::Vector4d(slice_plane.axis2.x(), slice_plane.axis2.y(), slice_plane.axis2.z(), 1.0)).block<3, 1>(0, 0);
		slice_plane.normal = slice_plane.axis1.cross(slice_plane.axis2);
		slice_plane.normal.normalize();
		Eigen::Vector3d c1 = slice_plane.axis1 + slice_plane.axis2;
		Eigen::Vector3d c2 = -slice_plane.axis1 + slice_plane.axis2;
		Eigen::Vector3d c3 = -slice_plane.axis1 - slice_plane.axis2;
		Eigen::Vector3d c4 =  slice_plane.axis1 - slice_plane.axis2;
		slice_plane.vertices.row(0) = slice_plane.center + 2.0 * c1;
		slice_plane.vertices.row(1) = slice_plane.center + 2.0 * c2;
		slice_plane.vertices.row(2) = slice_plane.center + 2.0 * c3;
		slice_plane.vertices.row(3) = slice_plane.center + 2.0 * c4;
		slice_plane.stored_vertices = slice_plane.vertices;

		viewer.slice_plane.set_vertices(slice_plane.vertices);
		updateVisibleIndices(slice_plane, V, T, tF, tC);
		viewer.erase_mesh(0);
		viewer.data().set_mesh(tV, tF);
		viewer.data().set_colors(tC);
		viewer.append_mesh();
		return true;
	}
	return false;
}


bool mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
	if (slice_plane.active) {
		down_rotation = trackball_angle;
		return true;
	}
	return false;
}


int main(int argc, char *argv[])
{
  // Load a mesh in MESH format
  igl::readMESH(TUTORIAL_SHARED_PATH "/cube.mesh", V, T);
  igl::convertMESH(V, T, tV, tF, tC);

  // Init the viewer
  igl::opengl::glfw::Viewer viewer;

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
  {
	  // Draw parent menu content
	  menu.draw_viewer_menu();

	  // Add slice plane
	  if (ImGui::CollapsingHeader("Slice Plane", ImGuiTreeNodeFlags_DefaultOpen))
	  {
		  static float dist = 0.f;
		  if (ImGui::Checkbox("enable", &slice_plane.enabled))
		  {
			  viewer.slice_enabled = slice_plane.enabled;
			  updateVisibleIndices(slice_plane, V, T, tF, tC);
			  viewer.erase_mesh(0);
			  viewer.data().set_mesh(tV, tF);
			  viewer.data().set_colors(tC);
			  viewer.append_mesh();
		  }

		  if (ImGui::Checkbox("active", &slice_plane.active))
		  {
			  updateVisibleIndices(slice_plane, V, T, tF, tC);
			  viewer.erase_mesh(0);
			  viewer.data().set_mesh(tV, tF);
			  viewer.data().set_colors(tC);
			  viewer.append_mesh();
		  }

		  static igl::opengl::glfw::Viewer::SlicePlaneOrientation dir = igl::opengl::glfw::Viewer::z;


		  if (ImGui::Combo("Orientation", (int *)(&dir), "x\0\y\0\z\0\0"))
		  {
			  switch (dir)
			  {
				  case igl::opengl::glfw::Viewer::x:
					  slice_plane.normal << 1.0, 0.0, 0.0;
					  slice_plane.axis1 << 0., 1., 0.;
					  slice_plane.axis2 << 0., 0., 1.;
					  slice_plane.vertices <<
						  0., -2., 2.,
						  0., 2., 2.,
						  0., 2., -2.,
						  0., -2., -2.;
					  slice_plane.stored_vertices = slice_plane.vertices;
					  break;
				  case igl::opengl::glfw::Viewer::y:
					  slice_plane.normal << 0.0, 1.0, 0.0;
					  slice_plane.axis1 << 0., 0., 1.;
					  slice_plane.axis2 << 1., 0., 0.;
					  slice_plane.vertices <<
						  -2., 0., 2.,
						  2., 0., 2.,
						  2., 0., -2.,
						  -2., 0., -2.;
					  slice_plane.stored_vertices = slice_plane.vertices;
					  break;
				  case igl::opengl::glfw::Viewer::z:
					  slice_plane.normal << 0.0, 0.0, 1.0;
					  slice_plane.axis1 << 1., 0., 0.;
					  slice_plane.axis2 << 0., 1., 0.;
					  slice_plane.vertices <<
						  -2., 2., 0.,
						  2, 2., 0.,
						  2., -2., 0.,
						  -2., -2., 0.;
					  slice_plane.stored_vertices = slice_plane.vertices;
					  break;
			  }
			  slice_plane.center = Eigen::Vector3d::Zero();
			  dist = 0.;
			  viewer.slice_plane.set_vertices(slice_plane.vertices);
			  updateVisibleIndices(slice_plane, V, T, tF, tC);
			  viewer.erase_mesh(0);
			  viewer.data().set_mesh(tV, tF);
			  viewer.data().set_colors(tC);
			  viewer.append_mesh();
		  }
		  

		  if (ImGui::SliderFloat("trans", &dist, -2.0f, 2.0f)) {
			  slice_plane.center = slice_plane.normal * (double)dist;
			  slice_plane.vertices = slice_plane.stored_vertices + (slice_plane.normal.transpose() * (double)dist).replicate(4, 1);
			  viewer.slice_plane.set_vertices(slice_plane.vertices);
			  updateVisibleIndices(slice_plane, V, T, tF, tC);
			  viewer.erase_mesh(0);
			  viewer.data().set_mesh(tV, tF);
			  viewer.data().set_colors(tC);
			  viewer.append_mesh();
		  }
	  }
  };

  // Register callbacks
  viewer.callback_mouse_move = &mouse_move;

  // Plot the mesh
  viewer.data().set_mesh(tV, tF);
  viewer.data().set_colors(tC);
  viewer.append_mesh();
  viewer.launch();
}

