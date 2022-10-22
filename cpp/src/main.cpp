#include <fstream>
#include <iostream>
#include <triangle/triangle.h>
#include <vector>
#include <Eigen/Dense>

extern "C" {
void exactinit();
uint32_t* triangulate_polygon(
    const double* points, const uint32_t* indices, const std::size_t* seperator, const std::size_t n_loops,
    const double* x_axis_data, const double* y_axis_data, const double* origin_data, uint32_t* n_triangles
);
}

static void writeOBJ(const std::string& name, const double* data, const uint32_t n_points, const std::vector<uint32_t>& indices) {
    std::ofstream file(name);
    for (uint32_t i = 0; i < n_points; i++) {
        file << "v " << data[2 * i] << " " << data[2 * i + 1] << " 0\n";
    }
    const auto n_triangles = indices.size() / 3;
    for (uint32_t i = 0; i < n_triangles; i++) {
        file << "f " << indices[3 * static_cast<std::vector<uint32_t, std::allocator<uint32_t>>::size_type>(i)] + 1
             << " " << indices[3 * i + 1] + 1 << " " << indices[3 * i + 2] + 1 << "\n";
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
    // std::vector<double> point_data{
    //     0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.25, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75,
    // };
    // auto result = triangualte(&point_data[0], static_cast<uint32_t>(point_data.size() / 2), nullptr, 0);
    const Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> point_data = Eigen::MatrixXd::Random(30000, 2);
    auto result = triangualte(&point_data(0, 0), static_cast<uint32_t>(point_data.size() / 2), nullptr, 0);
    writeOBJ("123.obj", &point_data(0, 0), point_data.rows(), result);
}