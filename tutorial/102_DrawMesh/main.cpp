#include <igl/readMESH.h>
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
		  // Expose variable directly ...
		  static bool slice_enabled = false;
		  if (ImGui::Checkbox("enable", &slice_enabled))
		  {
			  viewer.slice_enabled = slice_enabled;
		  }

		  static igl::opengl::glfw::Viewer::SlicePlaneOrientation dir = igl::opengl::glfw::Viewer::z;
		  if (ImGui::Combo("Orientation", (int *)(&dir), "x\0\y\0\z\0\0"))
		  {
			  Eigen::MatrixXd slice_plane_vertices;
			  slice_plane_vertices.resize(4, 3);
			  switch (dir)
			  {
				  case igl::opengl::glfw::Viewer::x:
					  slice_plane_vertices <<
						  0., -10., 10.,
						  0., 10., 10.,
						  0., 10., -10.,
						  0., -10., -10.;
					  break;
				  case igl::opengl::glfw::Viewer::y:
					  slice_plane_vertices <<
						  -10., 0., 10.,
						  10., 0., 10.,
						  10., 0., -10.,
						  -10., 0., -10.;
					  break;
				  case igl::opengl::glfw::Viewer::z:
					  slice_plane_vertices <<
						  -10., 10., 0.,
						  10., 10., 0.,
						  10., -10., 0.,
						  -10., -10., 0.;
					  break;
			  }
			  viewer.slice_plane.set_vertices(slice_plane_vertices);
			  viewer.slice_plane.dirty = igl::opengl::MeshGL::DIRTY_POSITION;
		  }
	  }
  };

  // Plot the mesh
  viewer.data().set_mesh(tV, tF);
  viewer.data().set_colors(tC);
  viewer.append_mesh();
  viewer.launch();
}

