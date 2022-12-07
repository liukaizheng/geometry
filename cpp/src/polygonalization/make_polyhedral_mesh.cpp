#include <polygonalization/conforming_mesh.h>
#include <polygonalization/make_polyhedral_mesh.h>
#include <triangle/tetrahedron.h>


void make_polyhedral_mesh_from_triangles(
    const double* points, const uint32_t n_points, const uint32_t* triangles, const uint32_t n_triangles
) {
    Constraints constraints(triangles, n_triangles);
    const auto mesh = TetMesh::tetrahedralize(points, n_points, 1e-8);
    place_virtual_constraints(mesh, constraints);
    // clang-format off
    // 0: tet; 1: f0; 2: f1; 3: f2; 4: f3.
    std::array<std::vector<std::vector<uint32_t>>, 5> tet_map{{
        std::vector<std::vector<uint32_t>>(mesh.tets.size()),
        std::vector<std::vector<uint32_t>>(mesh.tets.size()),
        std::vector<std::vector<uint32_t>>(mesh.tets.size()),
        std::vector<std::vector<uint32_t>>(mesh.tets.size()),
        std::vector<std::vector<uint32_t>>(mesh.tets.size())
    }};
    // clang-format on
    
}
