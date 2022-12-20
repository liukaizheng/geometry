#include <polygonalization/bsp_complex.h>
#include <polygonalization/conforming_mesh.h>
#include <polygonalization/make_polyhedral_mesh.h>
#include <triangle/tetrahedron.h>


void make_polyhedral_mesh_from_triangles(
    const double* points, const uint32_t n_points, const uint32_t* triangles, const uint32_t n_triangles,
    std::vector<double>& out_points, std::vector<uint32_t>& out_faces, std::vector<uint32_t>& seperators
) {
    Constraints constraints(triangles, n_triangles);
    auto mesh = TetMesh::tetrahedralize(points, n_points, 1e-8);
    place_virtual_constraints(mesh, constraints);
    // clang-format off
    // 0: f0; 1: f1; 2: f2; 3: f3; 4: tet.
    std::array<std::vector<std::vector<uint32_t>>, 5> tet_map{{
        std::vector<std::vector<uint32_t>>(mesh.tets.size()),
        std::vector<std::vector<uint32_t>>(mesh.tets.size()),
        std::vector<std::vector<uint32_t>>(mesh.tets.size()),
        std::vector<std::vector<uint32_t>>(mesh.tets.size()),
        std::vector<std::vector<uint32_t>>(mesh.tets.size())
    }};
    // clang-format on
    insert_constraints(mesh, constraints, tet_map.data());
    BSPComplex complex{mesh, &constraints, std::move(tet_map)};

    for (uint32_t cid = 0; cid < complex.cells.size();) {
        if (complex.cells[cid].constraints.empty()) {
            cid += 1;
        } else {
            complex.split_cell(cid);
        }
    }

    complex.decide_color();
    complex.complex_partition();
    complex.extract_skin(out_points, out_faces, seperators);
}
