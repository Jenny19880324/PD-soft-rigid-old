#include <igl/readMESH.h>
#include <igl/trackball.h>
#include <igl/unproject.h>
#include <igl/unproject_ray.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/MeshGL.h>
#include <igl/opengl/MarqueeGL.h>
#include <igl/opengl/FanGL.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/serialize.h>
#include <igl/dqs.h>
#include <igl/rbc.h>
#include <igl/rotation_matrix_from_axis_and_angle.h>
#include <GLFW/glfw3.h>
#include "tutorial_shared_path.h"
#include <unordered_set>
#include <unordered_map>
#include <igl/fit_hinged_rigid_motion.h>
#include <igl/png/writePNG.h>


bool vertex_pick_enabled = false;
bool vertices_marquee_enabled = false;
bool vertices_marquee_deselect_enabled = false;
bool vertices_move_enabled = false;
bool vertices_rotate_enabled = false;
bool external_force_enabled = false;
bool floor_enabled = false;
bool gravity_enabled = false;
bool output_moving_constraints = false;
bool output_screenshot = false;

Eigen::VectorXi temp_b;
Eigen::MatrixXd temp_U;
Eigen::Matrix<double, Eigen::Dynamic, 3> temp_bc;
Eigen::Matrix<double, Eigen::Dynamic, 3> temp_bc_rotate_base;
Eigen::RowVector3d center_of_temp_bc;
Eigen::MatrixXi visible_F;
Eigen::MatrixXd visible_C;

Eigen::MatrixXd V, U, P;
Eigen::MatrixXi T;
Eigen::MatrixXi F;
Eigen::MatrixXd C;
Eigen::Matrix<double, Eigen::Dynamic, 3> bc;
Eigen::VectorXi b;
Eigen::VectorXi N, A;
std::vector<std::vector<int>> I;

int pressed_b;
int anim_f = 0;
double anim_t = 0.0;
double anim_t_dir = 0.033;
igl::RBCData rbc_data; 
igl::opengl::MarqueeGL marquee_gl;
igl::opengl::FanGL fan_gl;


struct Marquee
{
	bool one_corner_set  = false;
	bool two_corners_set = false;
	int from_x;
	int from_y;
	int to_x;
	int to_y;

} marquee;

struct Fan
{
	bool center_set = false;
	bool first_point_on_circumstance_set = false;
	bool second_point_on_circumstance_set = false;
	Eigen::Vector2f center;
	Eigen::Vector2f first_point_on_circumstance;
	Eigen::Vector2f second_point_on_circumstance;
} fan;



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

struct FloorPlane {
	bool enabled = false;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	FloorPlane() {
		V.resize(4, 3);
		V << -10., 0., 10.,
			10., 0., 10.,
			10., 0., -10.,
			-10., 0., -10.;

		F.resize(2, 3);
		F << 0, 1, 3,
			1, 2, 3;
	}
} floor_plane;

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
	Eigen::VectorXi &b)
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

		float t;
		if (viewer.core.orthographic) {
			// I don't understand why I need to times it with 2.
			t3 = d * 2;
			t4 << ratio * t3.x(), t3.y(), -1, 1.0;
		}
		else {
			t = dnear / abs(d.z());
			t3 = t * d;
			t4 << ratio * t3.x(), t3.y(), t3.z(), 1.0;
		}
		t4 = view_model_inv * t4;
		frustum_points[i * 2] << t4.x() / t4.w(), t4.y() / t4.w(), t4.z() / t4.w();


		if (viewer.core.orthographic) {
			t3 = d * 2;
			t4 << ratio * t3.x(), t3.y(), 1, 1.0;
		}
		else {
			t = dfar / abs(d.z());
			t3 = t * d;
			t4 << ratio * t3.x(), t3.y(), t3.z(), 1.0;
		}
		t4 = view_model_inv * t4;
		frustum_points[i * 2 + 1] << t4.x() / t4.w(), t4.y() / t4.w(), t4.z() / t4.w();
	}

	// frustum normals
	Eigen::Vector3f frustum_normals[4];
	for (int i = 0; i < 4; ++i) {
		frustum_normals[i] = (frustum_points[i * 2] - frustum_points[i * 2 + 1]).cross(
			frustum_points[(i * 2 + 2) % 8] - frustum_points[(i * 2 + 3) % 8]);
		if (viewer.core.orthographic) {
			frustum_normals[i] = (frustum_points[i * 2] - frustum_points[i * 2 + 1]).cross(
				frustum_points[(i * 2 + 2) % 8] - frustum_points[i * 2]);
		}
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
			b.conservativeResize(b.rows() + 1);
			b.tail(1) << i;
		}
	}
}

void selection_difference(
	const Eigen::MatrixX3d &minuend_bc,
	const Eigen::VectorXi &minuend,
	const Eigen::VectorXi &subtrahend,
	Eigen::VectorXi &difference_b,
	Eigen::MatrixX3d &difference_bc)
{
	std::unordered_map<int, Eigen::RowVector3d> bc_map;
	for (int i = 0; i < minuend.rows(); i++) {
		bc_map[minuend(i)] = minuend_bc.row(i);
	}
	// difference = minuend - subtrahend
	int *m = const_cast<int *>(minuend.data());
	int *s = const_cast<int *>(subtrahend.data());
	std::vector<int> v(minuend.rows() + subtrahend.rows());
	std::vector<int>::iterator it;

	std::sort(m, m + minuend.rows());
	std::sort(s, s + subtrahend.rows());

	it = std::set_difference(m, m + minuend.rows(), s, s + subtrahend.rows(), v.begin());
	v.resize(it - v.begin());

	difference_b.resize(v.size());
	difference_bc.resize(difference_b.rows(), Eigen::NoChange);
	for (int i = 0; i < v.size(); i++) {
		difference_b(i) = v[i];
		difference_bc.row(i) = bc_map[v[i]];
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
		slice_plane.vertices.row(0) = slice_plane.center + c1;
		slice_plane.vertices.row(1) = slice_plane.center + c2;
		slice_plane.vertices.row(2) = slice_plane.center + c3;
		slice_plane.vertices.row(3) = slice_plane.center + c4;
		slice_plane.stored_vertices = slice_plane.vertices;

		viewer.slice_plane.set_vertices(slice_plane.vertices);
		updateVisibleIndices(slice_plane, U, T, C, visible_F, visible_C);
		viewer.data().clear();
		viewer.data().set_mesh(U, visible_F);
		viewer.data().set_colors(visible_C);
		return true;
	}

	if (vertex_pick_enabled && temp_b.rows() > 0) {
		Eigen::Vector3d cm = Eigen::Vector3d::Zero();
		for (int i = 0; i < temp_b.rows(); ++i) {
			cm += U.row(temp_b(i));
		}
		cm /= temp_b.size();
		cm = U.row(pressed_b);
		Eigen::Vector3f normal = cm.cast<float>() - viewer.core.camera_eye;


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
			temp_bc.resize(temp_b.rows(), Eigen::NoChange);
			for (int i = 0; i < temp_b.size(); ++i) {
				U.row(temp_b(i)) += delta;
				temp_bc.row(i) = U.row(temp_b(i));
			}
			viewer.data().set_vertices(U);

			if (b.rows() > bc.rows()) {
				bc.conservativeResize(bc.rows() + temp_bc.rows(), 3);
			}
			viewer.data().move_points(temp_bc, Eigen::RowVector3d(1.0, 0.0, 0.0));
			bc.block(bc.rows() - temp_bc.rows(), 0, temp_bc.rows(), 3) = temp_bc;

			igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
			return true;
		}
	}
	
	if (vertices_rotate_enabled && temp_b.rows() > 0 && temp_bc.rows() > 0) {
		if (fan.center_set && (fan.first_point_on_circumstance_set == false)) {
			fan.first_point_on_circumstance(0) = mouse_x;
			fan.first_point_on_circumstance(1) = mouse_y;
			fan_gl.dirty = true;

			float v0_x = fan.center.x();
			float v0_y = fan.center.y();
			float v1_x = fan.first_point_on_circumstance.x();
			float v1_y = fan.first_point_on_circumstance.y();


			fan_gl.V3.resize(2, Eigen::NoChange);
			fan_gl.V3 << v0_x, v0_y,
				v1_x, v1_y;
			return true;
		}

		if (fan.center_set && fan.first_point_on_circumstance_set && (fan.second_point_on_circumstance_set == false))
		{
			fan.second_point_on_circumstance(0) = mouse_x;
			fan.second_point_on_circumstance(1) = mouse_y;
			fan_gl.dirty = true;

			float v0_x = fan.center.x();
			float v0_y = fan.center.y();
			float v1_x = fan.first_point_on_circumstance.x();
			float v1_y = fan.first_point_on_circumstance.y();
			float v2_x = fan.second_point_on_circumstance.x();
			float v2_y = fan.second_point_on_circumstance.y();

			fan_gl.V3.resize(3, Eigen::NoChange);
			fan_gl.V3 << v0_x, v0_y,
				v1_x, v1_y,
				v2_x, v2_y;

			Eigen::Vector2d op1 = (fan.first_point_on_circumstance - fan.center).cast<double>();
			Eigen::Vector2d op2 = (fan.second_point_on_circumstance - fan.center).cast<double>();
			op1.normalize(); op2.normalize();
			double angle = acos(op1.dot(op2));

			Eigen::Vector3d axis = center_of_temp_bc.transpose() - viewer.core.camera_eye.cast<double>();
			axis.normalize();

			Eigen::Matrix3d rot_mat = igl::rotation_matrix_from_axis_and_angle(axis, angle);
			for (int i = 0; i < temp_bc.rows(); i++) {
				temp_bc.row(i) = (rot_mat * temp_bc_rotate_base.row(i).transpose()).transpose();
				U.row(temp_b(i)) = temp_bc.row(i);
			}
			viewer.data().set_vertices(U);
			viewer.data().move_points(temp_bc, Eigen::RowVector3d(1.0, 0.0, 0.0));
			bc.block(bc.rows() - temp_bc.rows(), 0, temp_bc.rows(), 3) = temp_bc;
			igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
			return true;
		}
		return true;
	}
	
	if (vertices_move_enabled && (!marquee.one_corner_set) && temp_b.rows() > 0 && temp_bc.rows() > 0)
	{
		Eigen::Vector3d pressed_V = U.row(pressed_b);
		Eigen::Vector3f normal = pressed_V.cast<float>() - viewer.core.camera_eye;

		Eigen::Vector3f s, dir;
		Eigen::Vector3f world_pos;
		igl::unproject_ray(Eigen::Vector2f(mouse_x, viewer.core.viewport(3) - mouse_y), viewer.core.view, viewer.core.proj, viewer.core.viewport, s, dir);
		dir.normalize(); normal.normalize();
		Eigen::Vector4f world_pressed_V;
		world_pressed_V << pressed_V.cast<float>(), 1.0;
		world_pressed_V = viewer.core.model * world_pressed_V;
		if (intersect_plane(normal, world_pressed_V.head<3>(), viewer.core.camera_eye, dir, world_pos)) {
			Eigen::Vector4f w;
			w << world_pos, 1.0;
			w = viewer.core.model.inverse() * w;
			world_pos = w.head<3>();
			Eigen::Vector3d delta = world_pos.cast<double>() - pressed_V;
			temp_bc.resize(temp_b.rows(), Eigen::NoChange);
			for (int i = 0; i < temp_b.size(); ++i) {
				U.row(temp_b(i)) += delta;
				temp_bc.row(i) = U.row(temp_b(i));
			}

			viewer.data().set_vertices(U);
			viewer.data().move_points(temp_bc, Eigen::RowVector3d(1.0, 0.0, 0.0));
			bc.block(bc.rows() - temp_bc.rows(), 0, temp_bc.rows(), 3) = temp_bc;
			//igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
			return true;
		}
	}

	if (vertices_marquee_enabled && marquee.one_corner_set) {
		marquee.to_x = mouse_x;
		marquee.to_y = mouse_y;
		marquee.two_corners_set = true;
		marquee_gl.dirty = true;

		float w = viewer.core.viewport(2);
		float h = viewer.core.viewport(3);

		float from_x = 2. * marquee.from_x / w - 1.;
		float from_y = -2. * marquee.from_y / h + 1.;
		float to_x = 2. * marquee.to_x / w - 1.;
		float to_y = -2. * marquee.to_y / h + 1.;

		marquee_gl.V_vbo << from_x, from_y, 0.,
			to_x, from_y, 0.,
			to_x, to_y, 0.,
			from_x, to_y, 0.;

		viewer.data().remove_points(temp_bc);
		temp_b.resize(0);
		temp_bc.resize(0, 3);
		marquee_select_vertices(viewer, marquee, U, temp_b);
		temp_bc.resize(temp_b.size(), 3);
		for (int i = 0; i < temp_b.size(); i++) {
			temp_bc.row(i) = U.row(temp_b(i));
		}

		viewer.data().add_points(temp_bc, Eigen::RowVector3d(1.0, 0.0, 0.0));
		return true;
	}

	if (vertices_marquee_deselect_enabled && marquee.one_corner_set) {
		marquee.to_x = mouse_x;
		marquee.to_y = mouse_y;
		marquee.two_corners_set = true;
		marquee_gl.dirty = true;

		float w = viewer.core.viewport(2);
		float h = viewer.core.viewport(3);

		float from_x = 2. * marquee.from_x / w - 1.;
		float from_y = -2. * marquee.from_y / h + 1.;
		float to_x = 2. * marquee.to_x / w - 1.;
		float to_y = -2. * marquee.to_y / h + 1.;

		marquee_gl.V_vbo << from_x, from_y, 0.,
			to_x, from_y, 0.,
			to_x, to_y, 0.,
			from_x, to_y, 0.;

		Eigen::VectorXi deselect_b;
		marquee_select_vertices(viewer, marquee, U, deselect_b);

		Eigen::VectorXi remained_b;
		Eigen::MatrixX3d remained_bc;
		selection_difference(bc, b, deselect_b, remained_b, remained_bc);
		viewer.data().set_points(remained_bc, Eigen::RowVector3d(0.0, 1.0, 0.0));
		
		Eigen::VectorXi remained_temp_b;
		Eigen::MatrixX3d remained_temp_bc;
		selection_difference(temp_bc, temp_b, deselect_b, remained_temp_b, remained_temp_bc);
		viewer.data().add_points(remained_temp_bc, Eigen::RowVector3d(1.0, 0.0, 0.0));


		b = remained_b;
		bc = remained_bc;
		temp_b = remained_temp_b;
		temp_bc = remained_temp_bc;
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
		if (viewer.data().bc.size() > 1 && 
			anim_f < viewer.data().bc.size()) {

			Eigen::VectorXi b_ith_frame = viewer.data().b[anim_f];
			Eigen::MatrixX3d bc_ith_frame = viewer.data().bc[anim_f];

			bc.block(0, 0, bc_ith_frame.rows(), 3) << bc_ith_frame;
			// green means b is precomputed, it's not temp_b anymore.
			viewer.data().move_points(bc_ith_frame, Eigen::RowVector3d(0.0, 1.0, 0.0), 0);
		}
		igl::rbc_precomputation(V, U, T, N, V.cols(), b, rbc_data);
		igl::rbc_solve(bc, rbc_data, U);
		viewer.data().set_vertices(U);
		viewer.data().compute_normals();

		//screen shot
		if (output_screenshot) {
			// Allocate temporary buffers
			Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> pngR(1280, 800);
			Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> pngG(1280, 800);
			Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> pngB(1280, 800);
			Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> pngA(1280, 800);

			// Draw the scene in the buffers
			viewer.core.draw_buffer(
				viewer.data(), false, pngR, pngG, pngB, pngA);

			// Save it to a PNG
			std::string filename = TUTORIAL_SHARED_PATH + std::string("/../screenshot/") + std::to_string(anim_f) + ".png";
			std::cout << "filename = " << filename << std::endl;
			const char *fname = filename.c_str();
			igl::png::writePNG(pngR, pngG, pngB, pngA, filename);
		}
		anim_t += anim_t_dir;
		anim_f++;
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
	case 'H':
	case 'h':
		if ((temp_b.rows() > 0) && (temp_bc.rows() > 0)) {
			bc.conservativeResize(bc.rows() + temp_bc.rows(), Eigen::NoChange);
			b.conservativeResize(b.rows() + temp_b.rows());
			bc.block(bc.rows() - temp_bc.rows(), 0, temp_bc.rows(), 3) = temp_bc;
			b.tail(temp_b.rows()) << temp_b;

			// to save it in the scene
			int number_of_frames = viewer.data().bc.size();
			if (number_of_frames == 0) {
				// no saved constraints in the scene now.
				viewer.data().b.push_back(b);
				viewer.data().bc.push_back(bc);
			}
			else {
				// append the fixed point constraint 
				// to the constraints in the scene
				for (int i = 0; i < number_of_frames; i++)
				{
					viewer.data().b[i].conservativeResize(b.rows());
					viewer.data().bc[i].conservativeResize(bc.rows(), Eigen::NoChange);
					viewer.data().b[i].tail(temp_b.rows()) << temp_b;
					viewer.data().bc[i].block(bc.rows() - temp_bc.rows(), 0, temp_bc.rows(), 3) << temp_bc;
				}
			}

			viewer.data().remove_points(temp_bc);
			viewer.data().add_points(temp_bc, Eigen::RowVector3d(0., 1., 0.));
			igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
			temp_b.resize(0);
			temp_bc.resize(0, 3);
		}
		return true;
	}
	return false;
}

int main(int argc, char *argv[])
{
  // Load a mesh in MESH format
  igl::readMESH(TUTORIAL_SHARED_PATH "/cube.mesh", V, T, F, C, N, A);
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
			  viewer.data().set_points(bc, Eigen::RowVector3d(0., 1., 0.));

			  // Find the bounding box
			  Eigen::Vector3d m = V.colwise().minCoeff();
			  Eigen::Vector3d M = V.colwise().maxCoeff();

			  double x_scale = (M - m).x();
			  double y_scale = (M - m).y();
			  double z_scale = (M - m).z();

			  slice_plane.normal << 0.0, 0.0, 1.0;
			  slice_plane.axis1 << x_scale, 0., 0.;
			  slice_plane.axis2 << 0., y_scale, 0.;
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
			  // Find the bounding box
			  Eigen::Vector3d m = V.colwise().minCoeff();
			  Eigen::Vector3d M = V.colwise().maxCoeff();

			  double x_scale = (M - m).x();
			  double y_scale = (M - m).y();
			  double z_scale = (M - m).z();

			  switch (dir)
			  {
				  case igl::opengl::glfw::Viewer::x:
					  slice_plane.normal << 1.0, 0.0, 0.0;
					  slice_plane.axis1 << 0., y_scale, 0.;
					  slice_plane.axis2 << 0., 0., z_scale;
					  slice_plane.vertices <<
						  0., -y_scale, z_scale,
						  0., y_scale, z_scale,
						  0., y_scale, -z_scale,
						  0., -y_scale, -z_scale;
					  slice_plane.stored_vertices = slice_plane.vertices;
					  break;
				  case igl::opengl::glfw::Viewer::y:
					  slice_plane.normal << 0.0, 1.0, 0.0;
					  slice_plane.axis1 << 0., 0., z_scale;
					  slice_plane.axis2 << x_scale, 0., 0.;
					  slice_plane.vertices <<
						  -x_scale, 0., z_scale,
						  x_scale, 0., z_scale,
						  x_scale, 0., -z_scale,
						  -x_scale, 0., -z_scale;
					  slice_plane.stored_vertices = slice_plane.vertices;
					  break;
				  case igl::opengl::glfw::Viewer::z:
					  slice_plane.normal << 0.0, 0.0, 1.0;
					  slice_plane.axis1 << x_scale, 0., 0.;
					  slice_plane.axis2 << 0., y_scale, 0.;
					  slice_plane.vertices <<
						  -x_scale, y_scale, 0.,
						  x_scale, y_scale, 0.,
						  x_scale, -y_scale, 0.,
						  -x_scale, -y_scale, 0.;
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
		  

		  if (ImGui::SliderFloat("trans", &dist, -15.0f, 15.0f)) {
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
				marquee_gl.init();
				if (vertices_marquee_deselect_enabled) {
					vertices_marquee_deselect_enabled = false;
				}
				// should call marquee_gl.free() somewhere.
		  }
		  //ImGui::SameLine();
		  if (ImGui::Checkbox("marquee vertices deselect", &vertices_marquee_deselect_enabled)) {
			  marquee_gl.init();
			  if (vertices_marquee_deselect_enabled) {
				  vertices_marquee_enabled = false;
			  }
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
	temp_b.resize(0);
	temp_bc.resize(0, 3);
	igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
}
	  }


  };

  menu.callback_draw_custom_window = [&]()
  {
	  // Simulation panel
	// Define next window position + size
	  ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 0), ImGuiSetCond_FirstUseEver);
	  ImGui::SetNextWindowSize(ImVec2(200, 400), ImGuiSetCond_FirstUseEver);
	  ImGui::Begin(
		  "Simulation", nullptr
	  );

	  if (ImGui::Button("reset")) {
		  U = V;
		  viewer.data().set_mesh(U, F);
		  viewer.data().compute_normals();

		  if (viewer.data().bc.size() > 0) {
			  b = viewer.data().b[0];
			  bc = viewer.data().bc[0];
		  }
		  viewer.data().set_points(bc, Eigen::RowVector3d(0., 1., 0.));
		  //viewer.data().remove_points(bc);
		  //viewer.data().b.clear();
		  //viewer.data().bc.clear();
		  //b.resize(0, Eigen::NoChange);
		  //bc.resize(0, Eigen::NoChange);
		  //temp_b.resize(0, Eigen::NoChange);
		  //temp_bc.resize(0, Eigen::NoChange);
		  anim_f = 0;
		  anim_t = 0.;

		  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
	  }

	  ImGui::SameLine();
	  ImGui::Text("%d\n", anim_f);
	  
	  
	  // material panel
	  if (ImGui::Combo("Material", (int *)(&rbc_data.energy), "PD material\0\RBC\0"))
	  {
		  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
	  }

	  ImGui::PushItemWidth(-80);
	  if (ImGui::DragInt("max iterations", &rbc_data.max_iter, 1.0f, 1, 100, "%.0f")) {

	  }

	  if (ImGui::Checkbox("with dynamics", &rbc_data.with_dynamics)) {
		  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
	  }

	  if (ImGui::Checkbox("collision", &rbc_data.collision_enabled)) {

	  }

	  if (ImGui::Checkbox("floor", &floor_enabled)) {
		  static int floor_idx = 0;
		  if (floor_enabled) {

			  // Find the bounding box
			  Eigen::Vector3d m = V.colwise().minCoeff();
			  Eigen::Vector3d M = V.colwise().maxCoeff();

			  double x_scale = (M - m).x();
			  double y_scale = (M - m).y();
			  double z_scale = (M - m).z();


			  floor_plane.V << -x_scale, 0., z_scale,
				  x_scale, 0., z_scale,
				  x_scale, 0., -z_scale,
				  -x_scale, 0., -z_scale;

			  viewer.append_mesh();
			  floor_idx = viewer.data_list.size() - 1;
			  viewer.data().set_mesh(floor_plane.V, floor_plane.F);
			  viewer.selected_data_index--;
		  }
		  else {
			  viewer.erase_mesh(floor_idx);
		  }
	  }
	  ImGui::SameLine();
	  if (ImGui::DragFloat("floor y", &rbc_data.floor_y, 0.1, -10.0, 10.0)) {
		  double floor_y = rbc_data.floor_y;
		  viewer.selected_data_index++;

		  floor_plane.V(0, 1) = floor_y;
		  floor_plane.V(1, 1) = floor_y;
		  floor_plane.V(2, 1) = floor_y;
		  floor_plane.V(3, 1) = floor_y;
		  viewer.data().set_vertices(floor_plane.V);

		  viewer.selected_data_index--;
	  }

	  // Expose the same variable directly
	  if (ImGui::DragFloat("mu", &rbc_data.mu, 0.0, 0.0, 10.0)) {
		  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
	  }
	  ImGui::PopItemWidth();

	  if (ImGui::Checkbox("external force", &external_force_enabled)) {
	  }
	  if (ImGui::Checkbox("gravity", &gravity_enabled)) {
		  if (gravity_enabled) {
			  if (rbc_data.with_dynamics == false) {
				  rbc_data.with_dynamics = true;
				  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
			  }
		  }
		  else {
			  rbc_data.g = 0.;
		  }
	  }
	  ImGui::SameLine();
	  ImGui::PushItemWidth(-40);
	  if (ImGui::DragFloat("g", &rbc_data.g, 1.0, -1.0, 0)) {
		  if (gravity_enabled) {
		  }
		  else {
			  rbc_data.g = 0.;
		  }
	  }
	  ImGui::PopItemWidth();

	  ImGui::PushItemWidth(-120);
	  if (ImGui::Combo("Constraint", (int *)(&rbc_data.constraint), "hard\0\soft\0")) {
		  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
	  }
	  if (ImGui::DragFloat("soft constraint weight", &rbc_data.constraint_weight, 1.0, 0.0, 100.))
	  {
		  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
	  }
	  if (ImGui::DragFloat("collision weight", &rbc_data.collision_weight, 1.0, 0.0, 100.0))
	  {
	  }

	  if (ImGui::Combo("Bone Constraint", (int *)(&rbc_data.bone_constraint), "affine\0rigid\0constrained\0")) {

	  }
	  ImGui::PopItemWidth();
	  if (ImGui::Button("Output Moving Constraint")) {
		  // create b and bc by code because it's difficult to 
		  // specify by viewport operation
		  if (temp_b.rows() == 0 || temp_bc.rows() == 0) {
			  std::cerr << "select vertices first" << std::endl;
		  }
		  else {
			  Eigen::RowVector3d cm = Eigen::RowVector3d::Zero(); // center of selected vertices
			  for (int i = 0; i < temp_bc.rows(); i++)
			  {
				  cm += temp_bc.row(i);
			  }
			  cm /= temp_bc.rows();

			  int number_of_frames = 100;
			  Eigen::RowVector3d t;
			  Eigen::Matrix3d R;
			  for (int i = 0; i < number_of_frames; i++)
			  {
				  double alpha = M_PI / 2 * (double)i / (double)number_of_frames;
				  double beta = M_PI * (double)i / (double)number_of_frames;
				  t.x() = 12. * (cos(alpha) - 1.);
				  t.y() = 12. * sin(alpha);
				  t.z() = 0.;

				  R = igl::rotation_matrix_from_axis_and_angle(
					  Eigen::Vector3d(0., 0., 1.),
					  beta);
				  Eigen::MatrixXd bc_ith_frame(temp_bc.rows(), temp_bc.cols());
				  for (int j = 0; j < temp_bc.rows(); j++) {
					  Eigen::RowVector3d delta = temp_bc.row(j) - cm;
					  bc_ith_frame.row(j) = delta * R.transpose() + cm + t;
				  }
				  viewer.data().b.push_back(temp_b);
				  viewer.data().bc.push_back(bc_ith_frame);
			  }


			  Eigen::VectorXi b_first_frame = viewer.data().b[0];
			  Eigen::MatrixX3d bc_first_frame = viewer.data().bc[0];
			 /* Eigen::VectorXi b_first_frame = temp_b;
			  Eigen::MatrixX3d bc_first_frame = temp_bc + Eigen::RowVector3d(6., 0., 0).replicate(temp_bc.rows(), 1);
			  viewer.data().b.push_back(temp_b);
			  viewer.data().bc.push_back(bc_first_frame);*/

			  b.conservativeResize(b.rows() + b_first_frame.rows());
			  b.block(b.rows() - b_first_frame.rows(), 0, b_first_frame.rows(), 1) << b_first_frame;
			  bc.conservativeResize(bc.rows() + bc_first_frame.rows(), Eigen::NoChange);
			  bc.block(bc.rows() - bc_first_frame.rows(), 0, bc_first_frame.rows(), 3) << bc_first_frame;
			  temp_b.resize(0);
			  temp_bc.resize(0, Eigen::NoChange);
			  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
			  output_moving_constraints = true;
		  }
	  }
	  ImGui::End();

	  // Mesh panel
	  ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 400), ImGuiSetCond_FirstUseEver);
	  ImGui::SetNextWindowSize(ImVec2(200, 80), ImGuiSetCond_FirstUseEver);
	  ImGui::Begin(
		  "Mesh", nullptr
	  );

	  ImGui::Text("#V %d\n", V.rows());
	  ImGui::Text("#T %d\n", T.rows());
	  ImGui::End();

	  // Output panel
	  ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 480), ImGuiSetCond_FirstUseEver);
	  ImGui::SetNextWindowSize(ImVec2(200, 80), ImGuiSetCond_FirstUseEver);
	  ImGui::Begin(
		  "Output", nullptr
	  );

	  if (ImGui::Checkbox("screenshot", &output_screenshot)) {

	  }
	  ImGui::End();

  };

  // Register callbacks
  viewer.callback_mouse_move = &mouse_move;
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.callback_mouse_down = [](igl::opengl::glfw::Viewer &viewer, int button, int modifier)->bool
  {
	  if (vertex_pick_enabled) {
		  int vid = -1;
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
			  temp_b.resize(temp_b.rows() + 1);
			  temp_b.tail<1>() << vid;
			  b.conservativeResize(b.rows() + temp_b.rows());
			  b.tail(temp_b.rows()) << temp_b;
			  pressed_b = vid;
			  return true;

		  }
	  }

	  if (vertices_rotate_enabled) {
		  if (fan.center_set && 
			  (!fan.first_point_on_circumstance_set) && 
			  (!fan.second_point_on_circumstance_set) &&
			  modifier == GLFW_MOD_ALT) {
			  std::cout << "first point on circumstance" << std::endl;
			  fan.first_point_on_circumstance(0) = viewer.current_mouse_x;
			  fan.first_point_on_circumstance(1) = viewer.current_mouse_y;
			  fan.first_point_on_circumstance_set = true;
			  return true;
		  }

		  if (fan.center_set &&
			  fan.first_point_on_circumstance_set &&
			  (!fan.second_point_on_circumstance_set) &&
			  modifier == GLFW_MOD_ALT) {
			  std::cout << "second point on circumstance" << std::endl;
			  fan.second_point_on_circumstance(0) = viewer.current_mouse_x;
			  fan.second_point_on_circumstance(1) = viewer.current_mouse_y;
			  fan.second_point_on_circumstance_set = true;
			  return true;
		  }
	  }

	  if (vertices_marquee_enabled) {
		  if (temp_b.rows() > 0 && temp_bc.rows() > 0) {
			  if (modifier == GLFW_MOD_ALT) { // rotate the selected vertices
				  std::cout << "rotate_enabled" << std::endl;
				  vertices_rotate_enabled = true;
				  b.conservativeResize(b.rows() + temp_b.rows());
				  b.tail(temp_b.rows()) << temp_b;
				  bc.conservativeResize(bc.rows() + temp_bc.rows(), Eigen::NoChange);
				  bc.block(bc.rows() - temp_bc.rows(), 0, temp_bc.rows(), 3) << temp_bc;
				  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);

				  center_of_temp_bc = Eigen::RowVector3d::Zero();
				  temp_bc_rotate_base = temp_bc;
				  for (int i = 0; i < temp_bc.rows(); i++) {
					  center_of_temp_bc += temp_bc.row(i);
				  }
				  center_of_temp_bc /= temp_bc.rows();

				  fan_gl.init();
				  fan.center(0) = viewer.current_mouse_x;
				  fan.center(1) = viewer.current_mouse_y;

				  float w = viewer.core.viewport(2);
				  float h = viewer.core.viewport(3);
				  fan_gl.proj << 2.f / w, 0.f, -1.f,
					  0.f, -2.f / h, 1.f,
					  0.f, 0.f, 1.f;
				  return true;
			  }
			  // If click on one of the vertices selected, move the selected vertices
			  int vid = -1;
			  double x = viewer.current_mouse_x;
			  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
			  if (igl::unproject_onto_mesh(
				  Eigen::Vector2f(x, y),
				  viewer.core.view * viewer.core.model,
				  viewer.core.proj,
				  viewer.core.viewport,
				  U, F, vid)) {
				  for (int i = 0; i < temp_b.rows(); i++) {
					  if (vid == temp_b(i)) {
						  vertices_move_enabled = true;
						  pressed_b = vid;
						  b.conservativeResize(b.rows() + temp_b.rows());
						  b.tail(temp_b.rows()) << temp_b;
						  bc.conservativeResize(bc.rows() + temp_bc.rows(), Eigen::NoChange);
						  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
						  return true;
					  }
				  }
			  }
		  }
		  vertices_move_enabled = false;
		  marquee.one_corner_set = true;
		  marquee.from_x = viewer.current_mouse_x;
		  marquee.from_y = viewer.current_mouse_y;
		  return true;
	  }

	  if (vertices_marquee_deselect_enabled) {
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
		  b.conservativeResize(b.rows() - temp_b.rows());
		  bc.conservativeResize(bc.rows() - temp_bc.rows(), Eigen::NoChange);
		  viewer.data().remove_points(temp_bc);
		  temp_b.resize(0);
		  temp_bc.resize(0, 3);
		  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
		  return true;
	  }

	  if (vertices_rotate_enabled) {
		  if ((!fan.center_set) && (!fan.first_point_on_circumstance_set) && (!fan.second_point_on_circumstance_set))
		  {
			  std::cout << "fan center set" << std::endl;
			  fan.center_set = true;
		  }
		  if (fan.center_set && fan.first_point_on_circumstance_set) {
			  std::cout << "restart the fan" << std::endl;
			  fan.center_set = false;
			  fan.first_point_on_circumstance_set = false;
			  fan.second_point_on_circumstance_set = false;
			  fan_gl.V3.resize(0, Eigen::NoChange);
			  fan_gl.V_vbo.resize(0, Eigen::NoChange);
			  vertices_rotate_enabled = false;
			  return true;
		  }
	  }

	  if (vertices_marquee_enabled) {
		  marquee.one_corner_set = false;
		  marquee.two_corners_set = false;

		  if (vertices_move_enabled) {
			  b.conservativeResize(b.rows() - temp_b.rows());
			  bc.conservativeResize(bc.rows() - temp_bc.rows(), 3);
			  viewer.data().remove_points(temp_bc);
			  temp_b.resize(0);
			  temp_bc.resize(0, 3);
			  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
			  vertices_move_enabled = false;
		  }

		  return true;
	  }

	  if (vertices_marquee_deselect_enabled) {
		  marquee.one_corner_set = false;
		  marquee.two_corners_set = false;
		  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);
		  return true;
	  }

	  return false;
  };

  viewer.callback_post_draw = [](igl::opengl::glfw::Viewer &viewer)
  {
	  if ((vertices_marquee_enabled||vertices_marquee_deselect_enabled)
		  && marquee.two_corners_set) {
		  marquee_gl.bind();
		  marquee_gl.draw();
	  }

	  if (vertices_rotate_enabled && fan.center_set) {
		  fan_gl.bind();
		  fan_gl.draw();
	  }
	  return false;
  };

  // Precomputation
  //rbc_data.max_iter = 100;
  //rbc_data.with_dynamics = true;
  //rbc_data.h = 0.5;
  igl::rbc_precomputation(V, T, N, V.cols(), b, rbc_data);

  // Plot the mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.core.is_animating = false;
  viewer.core.animation_max_fps = 30.;
  viewer.launch();
}

