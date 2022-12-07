#include <polygonalization/conforming_mesh.h>
#include <polygonalization/make_polyhedral_mesh.h>
#include <triangle/tetrahedron.h>


void make_polyhedral_mesh_from_triangles(const double* points, const uint32_t n_points, const uint32_t* triangles, const uint32_t n_triangles) {
    Constraints constraints(triangles, n_triangles);
    const auto mesh = TetMesh::tetrahedralize(points, n_points, 1e-8);
    place_virtual_constraints(mesh, constraints);
}
