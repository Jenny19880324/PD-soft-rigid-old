#include <igl/readMESH.h>
#include <igl/trackball.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/MeshGL.h>
#include <igl/unproject_onto_mesh.h>
#include "tutorial_shared_path.h"

bool vertex_pick_enabled = false;

Eigen::MatrixXd V;
Eigen::MatrixXi T;
Eigen::MatrixXi F;
Eigen::MatrixXd C;
Eigen::MatrixXi visible_F;
Eigen::MatrixXd visible_C;


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
void updateTriVisibility(
	const SlicePlane &slice_plane,
	const Eigen::PlainObjectBase<DerivedV> &V,
	const Eigen::PlainObjectBase<DerivedT> &F,
	std::vector<bool> &tri_visibility)
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
	tri_visibility.resize(F.rows());
	for (int tri_i = 0; tri_i < F.rows(); ++tri_i) {
		tri_visibility[tri_i] =
			vert_visibility[F(tri_i, 0)] &&
			vert_visibility[F(tri_i, 1)] &&
			vert_visibility[F(tri_i, 2)];
	}
}


template <typename DerivedV, typename DerivedT, typename DerivedF, typename DerivedC>
void updateVisibleIndices(
	const SlicePlane &slice_plane,
	const Eigen::PlainObjectBase<DerivedV> &V,
	const Eigen::PlainObjectBase<DerivedT> &F,
	const Eigen::PlainObjectBase<DerivedC> &C,
	Eigen::PlainObjectBase<DerivedF> &visible_F,
	Eigen::PlainObjectBase<DerivedC> &visible_C)
{
	visible_F.resize(F.rows(), 3); // max out the size for now
	visible_C.resize(F.rows(), 4);
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

	std::vector<bool> tri_visibility;
	updateTriVisibility(slice_plane, V, F, tri_visibility);

	size_t tri_visible_idx = 0;

	for (size_t tri_i = 0; tri_i < F.rows(); ++tri_i) {
		if (!tri_visibility[tri_i]) {
			continue;
		}

		visible_F(tri_visible_idx + 0, 0) = F(tri_i, 0);
		visible_F(tri_visible_idx + 0, 1) = F(tri_i, 1);
		visible_F(tri_visible_idx + 0, 2) = F(tri_i, 2);

		visible_C.row(tri_visible_idx) = C.row(tri_i);
		++tri_visible_idx;
	}
	visible_F.conservativeResize(tri_visible_idx, Eigen::NoChange);
	visible_C.conservativeResize(tri_visible_idx, Eigen::NoChange);
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
		updateVisibleIndices(slice_plane, V, F, C, visible_F, visible_C);
		viewer.data().clear();
		viewer.data().set_mesh(V, visible_F);
		viewer.data().set_colors(visible_C);
		return true;
	}
	return false;
}


int main(int argc, char *argv[])
{
  // Load a mesh in MESH format
  igl::readMESH(TUTORIAL_SHARED_PATH "/cube.mesh", V, T, F, C);

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

	  // Add slice plane panel
	  if (ImGui::CollapsingHeader("Slice Plane", ImGuiTreeNodeFlags_DefaultOpen))
	  {
		  static float dist = 0.f;
		  if (ImGui::Checkbox("enable", &slice_plane.enabled))
		  {
			  viewer.slice_enabled = slice_plane.enabled;
			  updateVisibleIndices(slice_plane, V, F, C, visible_F, visible_C);
			  viewer.data().clear();
			  viewer.data().set_mesh(V, visible_F);
			  viewer.data().set_colors(visible_C);
		  }

		  if (ImGui::Checkbox("active", &slice_plane.active))
		  {
			  updateVisibleIndices(slice_plane, V, F, C, visible_F, visible_C);
			  viewer.data().clear();
			  viewer.data().set_mesh(V, visible_F);
			  viewer.data().set_colors(visible_C);
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
			  updateVisibleIndices(slice_plane, V, F, C, visible_F, visible_C);
			  viewer.data().clear();
			  viewer.data().set_mesh(V, visible_F);
			  viewer.data().set_colors(visible_C);
		  }
		  

		  if (ImGui::SliderFloat("trans", &dist, -2.0f, 2.0f)) {
			  slice_plane.center = slice_plane.normal * (double)dist;
			  slice_plane.vertices = slice_plane.stored_vertices + (slice_plane.normal.transpose() * (double)dist).replicate(4, 1);
			  viewer.slice_plane.set_vertices(slice_plane.vertices);
			  updateVisibleIndices(slice_plane, V, F, C, visible_F, visible_C);
			  viewer.data().clear();
			  viewer.data().set_mesh(V, visible_F);
			  viewer.data().set_colors(visible_C);
		  }
	  }
	
	  // Add pick panel
	  if (ImGui::CollapsingHeader("pick", ImGuiTreeNodeFlags_DefaultOpen)) {
		  if (ImGui::Checkbox("vertex", &vertex_pick_enabled)) {
		  }
	  }
  };

  // Register callbacks
  viewer.callback_mouse_move = &mouse_move;
  viewer.callback_mouse_down = [](igl::opengl::glfw::Viewer &viewer, int, int)->bool
  {
	  if (vertex_pick_enabled) {
		  int vid = -1;
		  Eigen::Vector3f bc;
		  // Cast a ray in the view direction starting from the mouse position
		  double x = viewer.current_mouse_x;
		  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
		  if (igl::unproject_onto_mesh(
			  Eigen::Vector2f(x, y),
			  viewer.core.view * viewer.core.model,
			  viewer.core.proj,
			  viewer.core.viewport,
			  V, F, vid)) {
			  std::cout << "picked vid = " << vid << std::endl;
			  return true;

		  }
	  }
	  return false;
  };

  // Plot the mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.launch();
}

