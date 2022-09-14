#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include <iostream>
using namespace std;
using namespace geometrycentral::surface;

int main(int argc, char **argv) {
  cout << "hello world!" << endl;

  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readSurfaceMesh("arm.obj");
  return 0;
}
