#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <polygonalization/make_polyhedral_mesh.h>
#include <predicates/generic_point.h>
#include <predicates/interval_number.h>
#include <predicates/predicates.h>
#include <triangle/tetrahedron.h>
#include <triangle/triangle.h>
#include <vector>

extern "C" {
uint32_t* triangulate_polygon_soup(
    const double* points, const uint32_t* edge_data, const double* axis_data, const uint32_t* seperator,
    const uint32_t n_polygon, uint32_t* n_triangles
);
}

static void writeOBJ(
    const std::string& name, const double* data, const uint32_t n_points, const uint32_t* indices,
    const uint32_t n_triangles
) {
    std::ofstream file(name);
    for (uint32_t i = 0; i < n_points; i++) {
        file << "v " << data[3 * i] << " " << data[3 * i + 1] << " " << data[3 * i + 2] << "\n";
    }
    for (uint32_t i = 0; i < n_triangles; i++) {
        file << "f " << indices[3 * static_cast<std::vector<uint32_t, std::allocator<uint32_t>>::size_type>(i)] + 1
             << " " << indices[3 * i + 1] + 1 << " " << indices[3 * i + 2] + 1 << "\n";
    }
    file.close();
}

static void write_xyz(const std::string& name, const double* data, const uint32_t n_points) {
    std::ofstream file(name);
    for (uint32_t i = 0; i < n_points; i++) {
        const double* p = &data[i << 1];
        file << p[0] << " " << p[1] << " 0\n";
    }
    file.close();
}

using SPMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using T = Eigen::Triplet<double>;



int main() {
    std::vector<T> triples{{1, 0, 22.0}, {2, 0, 7.0}, {0, 1, 3.0},  {2, 1, 5.0},
                           {4, 2, 14.0}, {2, 3, 2.0}, {1, 4, 17.0}, {4, 4, 8.0}};
    const uint32_t n = 5;
    SPMat m(n, n);
    m.setFromTriplets(triples.begin(), triples.end());
    m.coeffRef(0, 1) = 0.0;
    m.prune(0.0);
    for (int k = 0; k < m.outerSize(); ++k) {
        std::cout << "k: " << k << "\n";
        for (SPMat::InnerIterator it(m, k); it; ++it) {
            it.value();
            it.row();   // row index
            it.col();   // col index (here it is equal to k)
            it.index(); // inner index, here it is equal to it.row()
            std::cout << "(" << it.row() << ", " << it.col() << "): " << it.value() << "\n";
        }
    }
    return 0;
}

/*int main() {
    exactinit();
    // const uint32_t n_points = 10000;
    // auto points = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>::Random(n_points, 3).eval();
    // const auto tets = TetMesh::tetrahedralize(points.data(), n_points, 1e-6);
    const uint32_t n_points = 8;
    std::vector<double> points {
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        0, 0, 1,
        1, 0, 1,
        1, 1, 1,
        0, 1, 1,
    };
    std::vector<uint32_t> indices{0, 3, 2, 0, 4, 3, 3, 4, 7, 3, 7, 6, 3, 6, 2, 6, 5, 2, 2, 5, 1, 1, 4, 0, 1, 5, 4, 4, 5,
6, 4, 6, 7}; make_polyhedral_mesh_from_triangles(points.data(), n_points, indices.data(), 11); return 0;
}*/
