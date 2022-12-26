#include <triangle/triangle.h>
#include <triangle/triangulate_polygon.h>

#include <Eigen/Dense>
#include <cmath>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <vector>

template <typename Index>
std::vector<Index> unique_indices(const Index* indices, const Index len, std::vector<Index>& new_to_origin_map) {
    new_to_origin_map.clear();
    Index count = 0;
    std::unordered_map<Index, Index> map;
    map.reserve(len);
    std::vector<Index> result(len);
    std::transform(indices, indices + len, result.begin(), [&count, &map](const Index old) {
        auto iter = map.find(old);
        if (iter != map.end()) {
            return iter->second;
        } else {
            map.emplace(old, count++);
            return count - 1;
        }
    });
    new_to_origin_map.resize(count);
    for (auto&& pair : std::move(map)) {
        new_to_origin_map[pair.second] = pair.first;
    }
    return result;
}

#ifdef EMSCRIPTEN
extern "C" {
#endif
uint32_t* triangulate_polygon(
    const double* points, const uint32_t* segments, const uint32_t n_segments, const double* origin_data,
    const double* x_axis_data, const double* y_axis_data, uint32_t* n_triangles
) {
    std::vector<uint32_t> new_to_origin_map;
    const auto new_segment = unique_indices(segments, n_segments << 1, new_to_origin_map);
    using Vector3d = Eigen::Map<const Eigen::Vector3d>;
    const Vector3d x_axis{x_axis_data};
    const Vector3d y_axis{y_axis_data};
    const Vector3d origin{origin_data};
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;
    Matrix points_2d(new_to_origin_map.size(), 2);
    // convert 3d points to 2d points
    for (uint32_t i = 0; i < new_to_origin_map.size(); i++) {
        const Vector3d p{points + new_to_origin_map[i] * 3};
        const auto v = (p - origin).eval();
        double* p_2d = &points_2d(i, 0);
        p_2d[0] = v.dot(x_axis);
        p_2d[1] = v.dot(y_axis);
    }
    auto triangle_data = triangulate(
        &points_2d(0, 0), static_cast<uint32_t>(new_to_origin_map.size()), &new_segment[0], n_segments, n_triangles
    );
    for (uint32_t i = 0; i < *n_triangles; i++) {
        uint32_t* data = &triangle_data[i * 3];
        data[0] = new_to_origin_map[data[0]];
        data[1] = new_to_origin_map[data[1]];
        data[2] = new_to_origin_map[data[2]];
    }
    return triangle_data;
}

uint32_t* triangulate_polygon_soup(
    const double* points, const uint32_t* edge_data, const double* axis_data, const uint32_t* seperator,
    const uint32_t n_polygon, uint32_t* n_triangles
) {
    std::vector<uint32_t> triangle_data;
    const uint32_t* segment_ptr = edge_data;
    const double* axis_ptr = axis_data;
    for (uint32_t i = 0; i < n_polygon; i++) {
        const uint32_t double_n_segments = seperator[i + 1] - seperator[i];
        const double* origin = axis_ptr;
        axis_ptr += 3;
        const double* x_axis = axis_ptr;
        axis_ptr += 3;
        const double* y_axis = axis_ptr;
        uint32_t n_poly_tri = 0;
        const uint32_t* tri_data =
            triangulate_polygon(points, segment_ptr, double_n_segments >> 1, origin, x_axis, y_axis, &n_poly_tri);
        std::copy(tri_data, tri_data + n_poly_tri * 3, std::back_inserter(triangle_data));
        delete[] tri_data;
        segment_ptr += double_n_segments;
        axis_ptr += 3;
    }
    *n_triangles = static_cast<uint32_t>(triangle_data.size() / 3);
    auto result = std::make_unique<uint32_t[]>(triangle_data.size());
    std::copy(triangle_data.begin(), triangle_data.end(), &result[0]);
    return result.release();
}
#ifdef EMSCRIPTEN
}
#endif
