#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <triangle/triangle.h>
#include <vector>

extern "C" {
void exactinit();
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
    std::vector<double> point_data{
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        0, 0, 1,
        1, 0, 1,
        1, 1, 1,
        0, 1, 1,
    };
    std::vector<uint32_t> indices{
        3, 2, 2, 1, 1, 0, 0, 3,
        4, 5, 5, 6, 6, 7, 7, 4,
        0, 1, 1, 5, 5, 4, 4, 0,
        1, 2, 2, 6, 6, 5, 5, 1,
        3, 7, 7, 6, 6, 2, 2, 3,
        4, 7, 7, 3, 3, 0, 0, 4,
    };
    std::vector<uint32_t> seperator{0, 8, 16, 24, 32, 40, 48};
    std::vector<double> axes {
        0, 0, 0, 0, 1, 0, 1, 0, 0,
        0, 0, 1, 1, 0, 0, 0, 1, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 1, 0, 0, 0, 1,
        1, 1, 0, -1, 0,0, 0, 0, 1,
        0, 1, 0, 0, 0, 1, 0, 1, 0,
    };
    uint32_t n_triangles = 0;
    const auto result = triangulate_polygon_soup(
        point_data.data(), indices.data(), axes.data(), seperator.data(), static_cast<uint32_t>(seperator.size() - 1),
        &n_triangles
    );
    writeOBJ("123.obj", point_data.data(), point_data.size() / 3, result, n_triangles);
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