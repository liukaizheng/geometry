#include <fstream>
#include <iostream>
#include <triangle/triangle.h>
#include <vector>
#include <Eigen/Dense>
#include <numeric>

extern "C" {
void exactinit();
uint32_t* triangulate_polygon(
    const double* points, const uint32_t* indices, const std::size_t* seperator, const std::size_t n_loops,
    const double* x_axis_data, const double* y_axis_data, const double* origin_data, uint32_t* n_triangles
);
}

static void writeOBJ(const std::string& name, const double* data, const uint32_t n_points, const uint32_t* indices, const uint32_t n_triangles) {
    std::ofstream file(name);
    for (uint32_t i = 0; i < n_points; i++) {
        file << "v " << data[2 * i] << " " << data[2 * i + 1] << " 0\n";
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

// int main() {
//     std::vector<double> point_data{
//         0.0,  0.0,  0.0, 1.0,  0.0,  0.0, 1.0,  1.0,  0.0, 0.0,  1.0,  0.0,
//         0.25, 0.25, 0.0, 0.75, 0.25, 0.0, 0.75, 0.75, 0.0, 0.25, 0.75, 0.0,
//     };
//     std::vector<uint32_t> indices{0, 2, 1, 3, 7, 6, 5, 4};
//     std::vector<std::size_t> seperator{4, 8};
//     std::vector<double> x_axis{1.0, 0.0, 0.0};
//     std::vector<double> y_axis{0.0, 1.0, 0.0};
//     std::vector<double> origin{0.0, 0.0, 0.0};
//     uint32_t n_triangles = 0;
//     const auto result = triangulate_polygon(
//         point_data.data(), indices.data(), seperator.data(), 2, x_axis.data(), y_axis.data(), origin.data(),
//         &n_triangles
//     );
//     return 0;
// }
int main() {
    exactinit();
    uint32_t n_triangles = 0;
    // std::vector<double> point_data{
    //    0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 2.0, 0.0, 3.0, 0.0, 4, 0, 5, 0, 6, 0, 7, 0, 2, 0, 4, 0, 5, 0, 2, 0
    //};
    // auto result = triangulate(&point_data[0], static_cast<uint32_t>(point_data.size() / 2), nullptr, 0, &n_triangles);
    auto n_point = 4000000;

    auto point_data = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>::Random(n_point, 2).eval();
    /* std::vector<uint32_t> indices(n_point);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j) { const auto pi = &point_data(i, 0);
        const auto pj = &point_data(j, 0);
        return pi[0] < pj[0] || (pi[0] == pj[0] && pi[1] < pj[1]);
    });
    const auto it = std::unique(indices.begin(), indices.end(), [&](auto i, auto j) { const auto pi = &point_data(i, 0);
        const auto pj = &point_data(j, 0);
        return pi[0] == pj[0] && pi[1] == pj[1];
    });
    indices.erase(it, indices.end());
    n_point = indices.size();
    auto copy = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>(n_point, 2);
    for (unsigned i = 0; i < n_point; i++) {
        copy.row(i) = point_data.row(indices[i]);
    }
    point_data = copy;*/
    // std::cout << point_data;
    /* std::vector<Eigen::Vector3d> v_arr;
    ReadXYZ("456.xyz", v_arr);
    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> point_data(v_arr.size(), 2);
    for (Eigen::Index i = 0; i < point_data.rows(); i++) {
        point_data(i, 0) = v_arr[i][0];
        point_data(i, 1) = v_arr[i][1];
    }
    uint32_t start = 0;
    n_point = point_data.rows();*/

    auto result =
        triangulate(&point_data(0, 0), n_point/* static_cast<uint32_t>(point_data.size() / 2)*/, nullptr, 0, &n_triangles);
    // write_xyz("123.xyz", &point_data(0, 0), n_point);
    writeOBJ("123.obj", &point_data(0, 0), n_point, result, n_triangles);
}