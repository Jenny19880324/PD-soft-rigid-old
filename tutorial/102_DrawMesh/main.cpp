#include <igl/readMESH.h>
#include <igl/trackball.h>
#include <igl/unproject.h>
#include <igl/unproject_ray.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl2/draw_rectangular_marquee.h>
#include <imgui/imgui.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/MeshGL.h>
#include <igl/opengl/MarqueeGL.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/rbc.h>
#include "tutorial_shared_path.h"

bool vertex_pick_enabled = false;
bool vertices_marquee_enabled = false;

Eigen::MatrixXd V, Vf, Vb, U;
Eigen::MatrixXi T;
Eigen::MatrixXi F;
Eigen::MatrixXd C;
Eigen::MatrixXi visible_F;
Eigen::MatrixXd visible_C;
Eigen::Matrix<double, Eigen::Dynamic, 3> bc;
Eigen::VectorXi b;
double anim_t = 0.0;
double anim_t_dir = 0.033;
igl::RBCData rbc_data; 
igl::opengl::MarqueeGL marquee_gl;

struct Marquee
{
	bool one_corner_set  = false;
	bool two_corners_set = false;
	int from_x;
	int from_y;
	int to_x;
	int to_y;

} marquee;



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
	const Eigen::PlainObjectBase<DerivedC> &C,
	Eigen::PlainObjectBase<DerivedF> &visible_F,
	Eigen::PlainObjectBase<DerivedC> &visible_C)
{
	visible_F.resize(T.rows() * 4, 3); // max out the size for now
	visible_C.resize(T.rows() * 4, 4);
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
	updateTriVisibility(slice_plane, V, F, tet_visibility);

	size_t tet_visible_idx = 0;

	for (size_t tet_i = 0; tet_i < T.rows(); ++tet_i) {
		if (!tet_visibility[tet_i]) {
			continue;
		}

		visible_F(tet_visible_idx * 4 + 0, 0) = T(tet_i, 0); 
		visible_F(tet_visible_idx * 4 + 0, 1) = T(tet_i, 1);
		visible_F(tet_visible_idx * 4 + 0, 2) = T(tet_i, 3);

		visible_F(tet_visible_idx * 4 + 1, 0) = T(tet_i, 1);
		visible_F(tet_visible_idx * 4 + 1, 1) = T(tet_i, 2);
		visible_F(tet_visible_idx * 4 + 1, 2) = T(tet_i, 3);

		visible_F(tet_visible_idx * 4 + 2, 0) = T(tet_i, 2);
		visible_F(tet_visible_idx * 4 + 2, 1) = T(tet_i, 0);
		visible_F(tet_visible_idx * 4 + 2, 2) = T(tet_i, 3);

		visible_F(tet_visible_idx * 4 + 3, 0) = T(tet_i, 0);
		visible_F(tet_visible_idx * 4 + 3, 1) = T(tet_i, 2);
		visible_F(tet_visible_idx * 4 + 3, 2) = T(tet_i, 1);

		visible_C.row(tet_visible_idx * 4 + 0) = C.row(tet_i * 4 + 0);
		visible_C.row(tet_visible_idx * 4 + 1) = C.row(tet_i * 4 + 1);
		visible_C.row(tet_visible_idx * 4 + 2) = C.row(tet_i * 4 + 2);
		visible_C.row(tet_visible_idx * 4 + 3) = C.row(tet_i * 4 + 3);
		++tet_visible_idx;
	}
	visible_F.conservativeResize(tet_visible_idx * 4, Eigen::NoChange);
	visible_C.conservativeResize(tet_visible_idx * 4, Eigen::NoChange);
}


bool intersect_plane(
	const Eigen::Vector3f &normal,
	const Eigen::Vector3f &point_on_plane,
	const Eigen::Vector3f &ray_origin,
	const Eigen::Vector3f &ray_dir,
	Eigen::Vector3f &intersection
)
{
	if (abs(normal.dot(ray_dir)) > 1e-6) {
		float t = (point_on_plane - ray_origin).dot(normal) / ray_dir.dot(normal);
		intersection = ray_origin + ray_dir * t;
		return true;
	}
	return false;
}


template <typename DerivedV>
void marquee_select_vertices(
	igl::opengl::glfw::Viewer &viewer, 
	const Marquee &marquee, 
	const Eigen::PlainObjectBase<DerivedV> &V,
	std::vector<int> &selected_vertices)
{
	const float width = viewer.core.viewport(2);
	const float height = viewer.core.viewport(3);

	float from_x = (float)marquee.from_x;
	float from_y = (float)marquee.from_y;
	float to_x = (float)marquee.to_x;
	float to_y = (float)marquee.to_y;

	const Eigen::Matrix4f model = viewer.core.model;
	const Eigen::Matrix4f view = viewer.core.view;
	const Eigen::Matrix4f proj = viewer.core.proj;

	const Eigen::Matrix4f view_model = view * model;
	const Eigen::Matrix4f view_model_inv = view_model.inverse();
	float dnear = viewer.core.camera_dnear;
	float dfar = viewer.core.camera_dfar;
	const Eigen::Vector3f eye = viewer.core.camera_eye;
	float ratio = width / height;

	Eigen::Vector3f screen_points[4] = {
		Eigen::Vector3f(from_x, to_y,   1.),
		Eigen::Vector3f(from_x, from_y, 1.),
		Eigen::Vector3f(to_x,   from_y, 1.),
		Eigen::Vector3f(to_x,   to_y,   1.)
	};

	Eigen::Matrix3f mat; // screen space to canonical space
	mat << 2. / width, 0., -1.,
			0., -2. / height, 1.,
			0., 0., 1.;

	for (int i = 0; i < 4; ++i) {
		screen_points[i] = mat * screen_points[i];
		screen_points[i].z() = -1. / std::tan(viewer.core.camera_view_angle / 360.0f * M_PI);
	}

	// 8 points of frustum
	Eigen::Vector3f frustum_points[8];
	Eigen::Vector3f t3;
	Eigen::Vector4f t4;
	for (int i = 0; i < 4; i++) {
		Eigen::Vector3f d = screen_points[i];
		float t = dnear / abs(d.z());
		t3 = t * d;
		t4 << ratio * t3.x(), t3.y(), t3.z(), 1.0;
		t4 = view_model_inv * t4;
		frustum_points[i * 2] << t4.x() / t4.w(), t4.y() / t4.w(), t4.z() / t4.w();

		t = dfar / abs(d.z());
		t3 = t * d;
		t4 << ratio * t3.x(), t3.y(), t3.z(), 1.0;
		t4 = view_model_inv * t4;
		frustum_points[i * 2 + 1] << t4.x() / t4.w(), t4.y() / t4.w(), t4.z() / t4.w();
	}

	// frustum normals
	Eigen::Vector3f frustum_normals[4];
	for (int i = 0; i < 4; ++i) {
		frustum_normals[i] = (frustum_points[i * 2] - frustum_points[i * 2 + 1]).cross(
			frustum_points[(i * 2 + 2) % 8] - frustum_points[(i * 2 + 3) % 8]);
		frustum_normals[i].normalize();
	}

	// Check which points are within the frustum formed by the selected rectangle
	for (int i = 0; i < V.rows(); ++i) {
		Eigen::Vector3f p = V.row(i).cast<float>();
		int j;
		for (j = 0; j < 4; ++j) {
			float dot_product = (p - frustum_points[j * 2]).dot(frustum_normals[j]);
			if (dot_product >= 0) {
				break;
			}
		}

		if (j == 4) {
			selected_vertices.push_back(i);
		}
	}
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
		updateVisibleIndices(slice_plane, U, T, C, visible_F, visible_C);
		viewer.data().clear();
		viewer.data().set_mesh(U, visible_F);
		viewer.data().set_colors(visible_C);
		return true;
	}

	if (vertex_pick_enabled && b.rows() > 0) {
		Eigen::Vector3d cm = Eigen::Vector3d::Zero();
		for (int i = 0; i < b.rows(); ++i) {
			cm += U.row(b(i));
		}
		cm /= b.size();
		Eigen::Vector3f normal = cm.cast<float>() - viewer.core.camera_eye ;

		Eigen::Vector3f s, dir;
		Eigen::Vector3f world_pos;
		igl::unproject_ray(Eigen::Vector2f(mouse_x, viewer.core.viewport(3) - mouse_y), viewer.core.view, viewer.core.proj, viewer.core.viewport, s, dir);
		dir.normalize(); normal.normalize();
		Eigen::Vector4f world_cm;
		world_cm << cm.cast<float>(), 1.0;
		world_cm = viewer.core.model * world_cm;
		if (intersect_plane(normal, world_cm.head<3>(), viewer.core.camera_eye, dir, world_pos)) {
			Eigen::Vector4f w;
			w << world_pos, 1.0;
			w = viewer.core.model.inverse() * w;
			world_pos = w.head<3>();
			Eigen::Vector3d delta = world_pos.cast<double>() - cm;
			bc.resize(b.size(), 3);
			for (int i = 0; i < b.size(); ++i) {
				U.row(b(i)) += delta;
				bc.row(i) = U.row(b(i));
			}
			viewer.data().set_vertices(U);
			viewer.data().set_points(bc, Eigen::RowVector3d(1.0, 0.0, 0.0));
			igl::rbc_precomputation(V, Vb, T, V.cols(), b, rbc_data);
			return true;
		}
	}
	
	if (vertices_marquee_enabled && marquee.one_corner_set) {
		marquee.to_x = mouse_x;
		marquee.to_y = mouse_y;
		marquee.two_corners_set = true;
		marquee_gl.dirty = true;

		float w= viewer.core.viewport(2);
		float h = viewer.core.viewport(3);

		float from_x = 2. * marquee.from_x / w - 1.;
		float from_y = -2. * marquee.from_y / h + 1.;
		float to_x = 2. * marquee.to_x / w - 1.;
		float to_y = -2. * marquee.to_y / h + 1.;

		marquee_gl.V_vbo << from_x, from_y, 0.,
							to_x, from_y, 0.,
							to_x, to_y, 0.,
							from_x, to_y, 0.;

		std::vector<int> selected_vertices;
		marquee_select_vertices(viewer, marquee, U, selected_vertices);
		Eigen::MatrixXd selected_V(selected_vertices.size(), 3);
		for (int i = 0; i < selected_vertices.size(); i++) {
			selected_V.row(i) = V.row(selected_vertices[i]);
		}
		viewer.data().add_points(selected_V, Eigen::RowVector3d(1., 0., 0.));
		return true;
	}
	return false;
}


bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
	using namespace Eigen;
	using namespace std;

	if (viewer.core.is_animating)
	{
		igl::rbc_solve(bc, rbc_data, U);
		viewer.data().set_vertices(U);
		viewer.data().compute_normals();

		anim_t += anim_t_dir;
	}
	return false;
}


bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
	switch (key)
	{
	case ' ':
		viewer.core.is_animating = !viewer.core.is_animating;
		return true;
	}
	return false;
}


int main(int argc, char *argv[])
{
  // Load a mesh in MESH format
  igl::readMESH(TUTORIAL_SHARED_PATH "/cube.mesh", V, Vf, Vb, T, F, C);
  U = V;

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
			  updateVisibleIndices(slice_plane, U, T, C, visible_F, visible_C);
			  viewer.data().clear();
			  viewer.data().set_mesh(U, visible_F);
			  viewer.data().set_colors(visible_C);
		  }

		  if (ImGui::Checkbox("active", &slice_plane.active))
		  {
			  updateVisibleIndices(slice_plane, U, T, C, visible_F, visible_C);
			  viewer.data().clear();
			  viewer.data().set_mesh(U, visible_F);
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
			  updateVisibleIndices(slice_plane, U, T, C, visible_F, visible_C);
			  viewer.data().clear();
			  viewer.data().set_mesh(U, visible_F);
			  viewer.data().set_colors(visible_C);
		  }
		  

		  if (ImGui::SliderFloat("trans", &dist, -2.0f, 2.0f)) {
			  slice_plane.center = slice_plane.normal * (double)dist;
			  slice_plane.vertices = slice_plane.stored_vertices + (slice_plane.normal.transpose() * (double)dist).replicate(4, 1);
			  viewer.slice_plane.set_vertices(slice_plane.vertices);
			  updateVisibleIndices(slice_plane, U, T, C, visible_F, visible_C);
			  viewer.data().clear();
			  viewer.data().set_mesh(U, visible_F);
			  viewer.data().set_colors(visible_C);
		  }
	  }
	
	  // Add pick panel
	  if (ImGui::CollapsingHeader("pick", ImGuiTreeNodeFlags_DefaultOpen)) {
		  if (ImGui::Checkbox("single vertex", &vertex_pick_enabled)) {
		  }
		  if (ImGui::Checkbox("marquee vertices", &vertices_marquee_enabled)) {
			  static bool initialized = false;
			  if (!initialized) {
				  marquee_gl.init();
				  initialized = true;
			  }
				// should call marquee_gl.free() somewhere.
		  }

		  if (ImGui::Button("clear")) {
			  marquee.one_corner_set = false;
			  marquee.two_corners_set = false;
			  marquee.from_x = -1;
			  marquee.from_y = -1;
			  marquee.to_x = -1;
			  marquee.to_y = -1;
			  viewer.data().points.resize(0, 0);
			  b.resize(0);
			  bc.resize(0, 3);
			  igl::rbc_precomputation(V, Vb, T, V.cols(), b, rbc_data);
		  }
	  }

	  // material panel
	  rbc_data.energy = igl::RBC_ENERGY_TYPE_RBC;
	  if (ImGui::Combo("Material", (int *)(&rbc_data.energy), "PD material\0\RBC\0"))
	  {
		  igl::rbc_precomputation(V, Vb, T, V.cols(), b, rbc_data);
	  }
  };

  menu.callback_draw_custom_window = [&]()
  {
	// Define next window position + size
	  ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 0), ImGuiSetCond_FirstUseEver);
	  ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
	  ImGui::Begin(
		  "Simulation", nullptr
	  );

	  // Expose the same variable directly
	  ImGui::PushItemWidth(-80);
	  ImGui::DragFloat("mu", &rbc_data.mu, 0.0, 0.0, 3.0);
	  ImGui::PopItemWidth();

	  ImGui::End();
  };

  // Register callbacks
  viewer.callback_mouse_move = &mouse_move;
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
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
			  U, F, vid)) {
			  viewer.data().add_points(U.row(vid), Eigen::RowVector3d(1., 0., 0.));
			  b.resize(b.rows() + 1);
			  b.tail<1>() << vid;
			  return true;

		  }
	  }

	  if (vertices_marquee_enabled) {
		  marquee.one_corner_set = true;
		  marquee.from_x = viewer.current_mouse_x;
		  marquee.from_y = viewer.current_mouse_y;
		  return true;
	  }
	  return false;
  };

  viewer.callback_mouse_up = [](igl::opengl::glfw::Viewer &viewer, int, int)->bool
  {
	  if (vertex_pick_enabled) {
		  viewer.data().points.resize(0, 0);
		  b.resize(0);
		  bc.resize(0, 3);
		  igl::rbc_precomputation(V, Vb, T, V.cols(), b, rbc_data);
		  return true;
	  }

	  if (vertices_marquee_enabled) {
		  marquee.one_corner_set = false;
		  marquee.two_corners_set = false;
		  return true;
	  }
	  return false;
  };

  viewer.callback_post_draw = [](igl::opengl::glfw::Viewer &viewer)
  {
	  if (vertices_marquee_enabled && marquee.two_corners_set) {
		  marquee_gl.bind();
		  marquee_gl.draw();
	  }
	  return false;
  };

  // Precomputation
  //rbc_data.max_iter = 100;
  rbc_data.with_dynamics = true;
  //rbc_data.h = 0.033;
  igl::rbc_precomputation(V, Vb, T, V.cols(), b, rbc_data);

  // Plot the mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.core.is_animating = false;
  viewer.core.animation_max_fps = 30.;
  viewer.launch();
}

