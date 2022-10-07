#include <Eigen/Dense>
#include <cmath>
#include <ranges>
#include <triangle/triangle.h>
#include <unordered_map>
#include <vector>

template <typename Index>
std::vector<Index> unique_indices(const Index* indices, const std::size_t len, std::vector<Index>& new_to_origin_map) {
    new_to_origin_map.clear();
    Index count = 0;
    std::unordered_map<Index, Index> map;
    map.reserve(len);
    std::vector<Index> result(static_cast<std::size_t>(len));
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

extern "C" {
uint32_t* triangulate_polygon(
    const double* points, const uint32_t* index_data, const std::size_t* seperator, const std::size_t n_loops,
    const double* x_axis_data, const double* y_axis_data, const double* origin_data, uint32_t* n_triangles
) {
    std::unordered_map<uint32_t, uint32_t> map;
    std::vector<uint32_t> new_to_origin_map;
    const auto indices = unique_indices(index_data, seperator[n_loops - 1], new_to_origin_map);
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

    const auto cross = [](const double* v1, const double* v2) -> double { return v1[0] * v2[1] - v1[1] * v2[0]; };
    const auto outer_loop_is_positive = [&seperator, &points_2d, &indices, &cross]() -> bool {
        double area = 0;
        for (std::size_t i = 0, len = seperator[0]; i < len; i++) {
            area += cross(&points_2d(indices[i], 0), &points_2d(indices[(i + 1) % len], 0));
        }
        return area > 0;
    };
    const auto v1_and_v2 = [](const double* p1, const double* p2, const double* p3) -> auto{
        return std::make_pair(
            Eigen::Vector2d{p3[0] - p2[0], p3[1] - p2[1]}.normalized(),
            Eigen::Vector2d{p1[0] - p2[0], p1[1] - p2[1]}.normalized()
        );
    };
    const auto hole_point = [&outer_loop_is_positive, &points_2d, &indices, &v1_and_v2,
                             &cross](const std::size_t start_idx, const std::size_t end_idx, double* data) {
        const bool positive = outer_loop_is_positive();
        const std::size_t len = end_idx - start_idx;
        std::size_t index = len;
        double min_angle_dot = 1.0;
        const static double EPSILON = std::cos(M_PI / 180.0);
        for (std::size_t i = 0; i < len; i++) {
            const auto view{
                std::vector{len + i - 1, i, i + 1} |
                std::views::transform([len, &points_2d, &indices, start_idx](const auto idx) {
                    return &points_2d(indices[idx % len + start_idx], 0);
                })};
            Eigen::Vector2d v1, v2;
            std::tie(v1, v2) = v1_and_v2(view[0], view[1], view[2]);
            if (positive ^ (cross(&v1[0], &v2[0]) > 0)) {
                const double dot = std::abs(v1.dot(v2));
                if (dot < 1. - EPSILON) {
                    data[0] = (view[0][0] + view[2][0]) * 0.5;
                    data[1] = (view[0][1] + view[2][1]) * 0.5;
                    data[2] = (view[0][2] + view[2][2]) * 0.5;
                    return;
                } else {
                    if (dot < min_angle_dot) {
                        min_angle_dot = dot;
                        index = i;
                    }
                }
            }
        }
        const auto three{
            std::vector{len + index - 1, index, index + 1} |
            std::views::transform([len, &points_2d, &indices, start_idx](const auto idx) {
                return &points_2d(indices[idx % len + start_idx], 0);
            })};
        data[0] = (three[0][0] = three[2][0]) * 0.5;
        data[1] = (three[0][1] = three[2][1]) * 0.5;
        data[2] = (three[0][2] = three[2][2]) * 0.5;
    };
    std::vector<double> hole_points((n_loops - 1) * 2);
    for (std::size_t i = 1; i < n_loops; i++) {
        hole_point(seperator[i - 1], seperator[i], &hole_points[(i - 1) * 2]);
    }
    triangulateio in;

    // points
    in.numberofpoints = static_cast<int>(points_2d.rows());
    in.pointlist = &points_2d(0, 0);
    in.pointmarkerlist = nullptr;
    in.numberofpointattributes = 0;

    // triangles
    in.numberoftriangles = 0;
    in.trianglelist = nullptr;
    in.numberofcorners = 0;
    in.numberoftriangleattributes = 0;
    in.triangleattributelist = nullptr;

    // segments
    in.numberofsegments = static_cast<int>(seperator[n_loops - 1]);
    std::vector<int> edges;
    edges.reserve(seperator[n_loops - 1]);
    for (uint32_t idx = 0; idx < n_loops; idx++) {
        const std::size_t start_idx = idx == 0 ? 0 : seperator[idx - 1];
        const std::size_t end_idx = seperator[idx];
        const auto len = end_idx - start_idx;
        for (std::size_t i = 0; i < len; i++) {
            edges.emplace_back(static_cast<int>(indices[i + start_idx]));
            edges.emplace_back(static_cast<int>(indices[(i + 1) % len + start_idx]));
        }
    }
    in.segmentlist = edges.data();
    in.segmentmarkerlist = nullptr;

    // holes
    in.numberofholes = static_cast<int>(hole_points.size() >> 1);
    in.holelist = in.numberofholes == 0 ? nullptr : hole_points.data();
    in.numberofregions = 0;
    in.neighborlist = nullptr;

    triangulateio out;
    out.pointlist = nullptr;
    out.pointmarkerlist = nullptr;
    out.pointattributelist = nullptr;
    out.trianglelist = nullptr;
    out.segmentlist = nullptr;
    out.segmentmarkerlist = nullptr;
    out.edgemarkerlist = nullptr;
    out.holelist = nullptr;

    const std::string str{"QpzBPN"};
    triangulate(const_cast<char*>(str.data()), &in, &out, nullptr);

    *n_triangles = static_cast<uint32_t>(out.numberoftriangles);
    const auto n_triangle_arr = *n_triangles * 3;
    auto result = std::make_unique<uint32_t[]>(n_triangle_arr);
    for (uint32_t i = 0; i < n_triangle_arr; i++) {
        result[i] = new_to_origin_map[static_cast<std::size_t>(out.trianglelist[i])];
    }
    free(out.trianglelist);
    return result.release();
}
}
