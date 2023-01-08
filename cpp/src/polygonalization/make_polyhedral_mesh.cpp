#include <algorithm>
#include <numeric>
#include <unordered_map>

#include <polygonalization/bsp_complex.h>
#include <polygonalization/conforming_mesh.h>
#include <polygonalization/make_polyhedral_mesh.h>
#include <triangle/tetrahedron.h>
#include <triangle/triangulate_polygon.h>


void make_polyhedral_mesh_from_triangles(
    const double* points, const uint32_t n_points, const uint32_t* triangles, const uint32_t n_triangles, const std::vector<std::unordered_map<uint32_t, uint32_t>>& ori_edge_parents,
    std::vector<double>& out_points, std::vector<uint32_t>& out_faces, std::vector<double>& axes,
    std::vector<uint32_t>& seperators
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
    BSPComplex complex{mesh, &constraints, ori_edge_parents, std::move(tet_map)};

    for (uint32_t cid = 0; cid < complex.cells.size();) {
        if (complex.cells[cid].constraints.empty()) {
            cid += 1;
        } else {
            complex.split_cell(cid);
        }
    }

    complex.decide_color();
    complex.complex_partition();
    complex.extract_skin(out_points, out_faces, axes, seperators);
}

static void remove_duplicates(
    const double* points, const uint32_t n_points, std::vector<double>& out_points, std::vector<uint32_t>& pmap
) {
    constexpr double epsilon = 0.000001;
    int e;
    std::frexp(epsilon, &e);
    e -= 1;
    const double base = std::ldexp(1.0, e);
    std::vector<double> approx_points(n_points * 3);
    for (uint32_t i = 0; i < approx_points.size(); i++) {
        approx_points[i] = std::ldexp(std::round(points[i] / base), e);
    }
    std::vector<uint32_t> indices(n_points);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&approx_points](const uint32_t i, const uint32_t j) {
        const double* pi = &approx_points[i * 3];
        const double* pj = &approx_points[j * 3];
        if (pi[0] < pj[0]) {
            return true;
        } else if (pi[0] > pj[0]) {
            return false;
        } else {
            return pi[1] < pj[1] || (pi[1] == pj[1] && pi[2] < pj[2]);
        }
    });

    std::vector<uint32_t> out_pts_indices{indices[0]};
    pmap.resize(n_points, 0);
    for (uint32_t k = 1; k < n_points; k++) {
        const uint32_t i = indices[k];
        const uint32_t j = indices[k - 1];
        const double* pi = &approx_points[i * 3];
        const double* pj = &approx_points[j * 3];
        if (pi[0] == pj[0] && pi[1] == pj[1] && pi[2] == pj[2]) {
            pmap[i] = static_cast<uint32_t>(out_pts_indices.size() - 1);
        } else {
            pmap[i] = static_cast<uint32_t>(out_pts_indices.size());
            out_pts_indices.emplace_back(i);
        }
    }
    out_points.resize(out_pts_indices.size() * 3);
    for (uint32_t i = 0; i < out_pts_indices.size(); i++) {
        const double* src = &approx_points[out_pts_indices[i] * 3];
        double* tar = &out_points[i * 3];
        tar[0] = src[0];
        tar[1] = src[1];
        tar[2] = src[2];
    }
}

static decltype(auto) get_edge_parents(std::vector<uint32_t>& edge_data, const uint32_t n_points) {
    std::vector<std::unordered_map<uint32_t, uint32_t>> edge_parents(n_points - 1);
    for (uint32_t i = 0; i < edge_data.size(); i += 2) {
        uint32_t pi = edge_data[i];
        uint32_t pj = edge_data[i + 1];
        if (pi > pj) {
            std::swap(pi, pj);
        }
        auto& map = edge_parents[pi];
        auto it = map.find(pj);
        if (it == map.end()) {
            map.emplace(pj, i >> 1);
        }
    }
    return edge_parents;
}

extern "C" {
uint32_t make_polyhedral_mesh(
    const double* points, const uint32_t n_points, const uint32_t* edge_data, const double* axis_data,
    const uint32_t* seperator, const uint32_t n_polygons, double** out_points, uint32_t** out_polygons,
    double** out_axis_data, uint32_t** out_seperators
) {
    std::vector<double> unique_points;
    std::vector<uint32_t> pmap;
    remove_duplicates(points, n_points, unique_points, pmap);
    std::vector<uint32_t> edges(seperator[n_polygons]);
    for (uint32_t i = 0; i < edges.size(); i++) {
        edges[i] = pmap[edge_data[i]];
    }
    const auto edge_parents = get_edge_parents(edges, static_cast<uint32_t>(unique_points.size() - 1));
    
    uint32_t n_triangles;
    const uint32_t* triangles =
        triangulate_polygon_soup(unique_points.data(), edges.data(), axis_data, seperator, n_polygons, &n_triangles);
    std::vector<double> out_pts_vec;
    std::vector<uint32_t> out_polys_vec;
    std::vector<double> axes_vec;
    std::vector<uint32_t> out_seperator_vec;
    make_polyhedral_mesh_from_triangles(
        unique_points.data(), static_cast<uint32_t>(unique_points.size() / 3), triangles, n_triangles, edge_parents, out_pts_vec,
        out_polys_vec, axes_vec, out_seperator_vec
    );
    delete[] triangles;

    // points
    auto out_pts = std::make_unique<double[]>(out_pts_vec.size());
    std::copy(out_pts_vec.begin(), out_pts_vec.end(), out_pts.get());
    *out_points = out_pts.release();

    // polygons
    auto out_polys = std::make_unique<uint32_t[]>(out_polys_vec.size());
    std::copy(out_polys_vec.begin(), out_polys_vec.end(), out_polys.get());
    *out_polygons = out_polys.release();

    // axes
    auto axes = std::make_unique<double[]>(axes_vec.size());
    std::copy(axes_vec.begin(), axes_vec.end(), axes.get());
    *out_axis_data = axes.release();

    // seperators
    auto out_seps = std::make_unique<uint32_t[]>(out_seperator_vec.size());
    std::copy(out_seperator_vec.begin(), out_seperator_vec.end(), out_seps.get());
    *out_seperators = out_seps.release();

    return static_cast<uint32_t>(out_seperator_vec.size() - 1);
}
}
