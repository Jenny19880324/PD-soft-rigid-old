#include <igl/readMESH.h>
#include <igl/opengl/glfw/Viewer.h>
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

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(tV, tF);
  viewer.data().set_colors(tC);
  viewer.launch();
}

