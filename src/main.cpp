#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include <iostream>
using namespace std;
using namespace geometrycentral::surface;

int main(int argc, char **argv) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = loadMesh("cube.obj");
  for (int i = 0; i < 8; i++) {
    VertexData<bool> isOrigVert(*mesh, true);
    EdgeData<bool> isOrigEdge(*mesh, true);
    std::vector<Edge> toFlip;
    for (Edge edge : mesh->edges()) {
      if (!isOrigEdge[edge]) {
        continue;
      }
      const Vertex oldA = edge.halfedge().tipVertex();
      const Vertex oldB = edge.halfedge().tailVertex();
      const auto oldAPos = geometry->vertexPositions[oldA];
      const auto oldBPos = geometry->vertexPositions[oldB];

      const auto newV = mesh->splitEdgeTriangular(edge).vertex();
      const auto newPos = 0.5 * (oldAPos + oldBPos);
      geometry->vertexPositions[newV] = newPos;
      for (Edge e : newV.adjacentEdges()) {
        isOrigEdge[e] = false;
        Vertex otherV = e.otherVertex(newV);
        if (isOrigVert[otherV] && otherV != oldA && otherV != oldB) {
          toFlip.push_back(e);
        }
      }
    }
    for (const auto e : std::move(toFlip)) {
      mesh->flip(e);
    }
  }
  writeSurfaceMesh(*mesh, *geometry, "cube1.obj");
  return 0;
}
