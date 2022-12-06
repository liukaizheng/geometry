#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <triangle/triangle.h>
#include <triangle/tetrahedron.h>
#include <vector>
#include <predicates/interval_number.h>
#include <predicates/generic_point.h>
#include <predicates/predicates.h>

extern "C" {
uint32_t* triangulate_polygon_soup(
    const double* points, const uint32_t* edge_data, const double* axis_data, const uint32_t* seperator, const uint32_t n_polygon,
     uint32_t* n_triangles
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

int main() {
    exactinit();
    const uint32_t n_points = 10000;
    auto points = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>::Random(n_points, 3).eval();
    const auto tets = TetMesh::tetrahedralize(points.data(), n_points, 1e-6);
    // const uint32_t n_points = 8;
    // std::vector<double> points {
    //     0, 0, 0,
    //     1, 0, 0,
    //     1, 1, 0,
    //     0, 1, 0,
    //     0, 0, 1,
    //     1, 0, 1,
    //     1, 1, 1,
    //     0, 1, 1,
    // };
    // const auto tets = Tetrahedrons::tetrahedralize(points.data(), n_points, 1e-6);
    std::vector<uint32_t> indices;
    // for (const auto& tet : tets.tets) {
    //     for (const auto vid : {tet.data[0], tet.data[1], tet.data[2]}) {
    //         indices.emplace_back(vid);
    //     }
    //     for (const auto vid : {tet.data[0], tet.data[3], tet.data[1]}) {
    //         indices.emplace_back(vid);
    //     }
    //     for (const auto vid : {tet.data[1], tet.data[3], tet.data[2]}) {
    //         indices.emplace_back(vid);
    //     }
    //     for (const auto vid : {tet.data[2], tet.data[3], tet.data[0]}) {
    //         indices.emplace_back(vid);
    //     }
    // }
    for (const auto& tet : tets.tets) {
        if (tet.data[3] == tets.n_points) {
            indices.emplace_back(tet.data[0]);
            indices.emplace_back(tet.data[1]);
            indices.emplace_back(tet.data[2]);
        }
    }
    writeOBJ("123.obj", points.data(), n_points, indices.data(), indices.size() / 3);
    return 0;
}
// int main() {
//     exactinit();
//     uint32_t n_triangles = 0;

//     uint32_t n_point = 1000000;
//     auto point_data = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>::Random(n_point, 2).eval();
//     std::vector<uint32_t> indices(n_point);
//     std::iota(indices.begin(), indices.end(), 0);
//     std::sort(indices.begin(), indices.end(), [&](int i, int j) { const auto pi = &point_data(i, 0);
//         const auto pj = &point_data(j, 0);
//         return pi[0] < pj[0] || (pi[0] == pj[0] && pi[1] < pj[1]);
//     });
//     const auto it = std::unique(indices.begin(), indices.end(), [&](auto i, auto j) { const auto pi = &point_data(i, 0);
//         const auto pj = &point_data(j, 0);
//         return pi[0] == pj[0] && pi[1] == pj[1];
//     });
//     indices.erase(it, indices.end());
//     n_point = indices.size();
//     auto copy = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>(n_point, 2);
//     for (unsigned i = 0; i < n_point; i++) {
//         copy.row(i) = point_data.row(indices[i]);
//     }
//     point_data = copy;
//     std::vector<uint32_t> segment{0, n_point - 1};

//     auto start = std::chrono::high_resolution_clock::now();
//     auto result = triangulate(
//         &point_data(0, 0), static_cast<uint32_t>(point_data.size() / 2), &segment[0], segment.size() / 2, &n_triangles
//     );
//     std::cout << "time elapsed "
//               << std::chrono::duration_cast<std::chrono::milliseconds>(
//                      std::chrono::high_resolution_clock::now() - start
//                  )
//                      .count()
//               << " ms\n";
//     writeOBJ("123.obj", &point_data(0,0), n_point, result, n_triangles);
// }