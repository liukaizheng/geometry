#include <iostream>
#include <vector>

extern "C" {
uint32_t* triangulate_polygon(
    const double* points, const uint32_t* indices, const std::size_t* seperator, const std::size_t n_loops,
    const double* x_axis_data, const double* y_axis_data, const double* origin_data, uint32_t* n_triangles
);
}

int main() {
    std::vector<double> point_data{
        0.0,  0.0,  0.0, 1.0,  0.0,  0.0, 1.0,  1.0,  0.0, 0.0,  1.0,  0.0,
        0.25, 0.25, 0.0, 0.75, 0.25, 0.0, 0.75, 0.75, 0.0, 0.25, 0.75, 0.0,
    };
    std::vector<uint32_t> indices{0, 1, 2, 3, 7, 6, 5, 4};
    std::vector<std::size_t> seperator{4, 8};
    std::vector<double> x_axis{1.0, 0.0, 0.0};
    std::vector<double> y_axis{0.0, 1.0, 0.0};
    std::vector<double> origin{0.0, 0.0, 0.0};
    uint32_t n_triangles = 0;
    const auto result = triangulate_polygon(
        point_data.data(), indices.data(), seperator.data(), 2, x_axis.data(), y_axis.data(), origin.data(), &n_triangles
    );
    return 0;
}