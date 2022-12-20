#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <fstream>
#include <graphcut/graphcut.h>
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

static void read_flow(
    std::vector<double>& src_cap, std::vector<double>& sink_cap,
    std::vector<std::unordered_map<uint32_t, double>>& edges
) {
    uint32_t n_nodes = 498;
    uint32_t s = 1, t = 500;
    src_cap.resize(n_nodes, 0);
    sink_cap.resize(n_nodes, 0);
    edges.resize(n_nodes);
    std::ifstream in("max_flow9.dat");
    std::string str;
    auto index = [&](const uint32_t i) -> uint32_t {
        const uint32_t d = (i > s) + (i > t);
        return i - d - 1;
    };
    while (std::getline(in, str)) {
        std::stringstream ss(str);
        std::string _s;
        uint32_t r, c;
        double val;
        ss >> _s >> r >> c >> val;

        if (r == t || c == s) {
            continue;
        } else if (r == s) {
            src_cap[index(c)] = val;
        } else if (c == t) {
            sink_cap[index(r)] = val;
        } else {
            uint32_t a = index(r);
            uint32_t b = index(c);
            if (a > b) {
                std::swap(a, b);
            }
            auto& map = edges[a];
            if (map.find(b) == map.end()) {
                map.emplace(b, val);
            }
        }
    }
}

/* int main() {
    const uint32_t n_nodes = 1000;
    auto src = Eigen::VectorXd::Random(n_nodes).eval();
    auto sink = Eigen::VectorXd::Random(n_nodes).eval();
    auto cap = Eigen::MatrixXd::Random(n_nodes, n_nodes).eval();
    std::vector<double> srcs_cap(n_nodes, 0);
    for (uint32_t i = 0; i < 100; i++) {
        srcs_cap[i] = std::abs(std::round(src[i] * n_nodes));
    }
    std::vector<double> sink_cap(n_nodes, 0);
    for (uint32_t i = 0; i < 100; i++) {
         sink_cap[n_nodes - 1 - i] = std::abs(std::round(sink[i] * n_nodes));
    }
    GraphCut g(n_nodes, srcs_cap.data(), sink_cap.data());
    for (uint32_t i = 0; i < n_nodes; i++) {
        for (uint32_t j = i + 1; j < n_nodes; j++) {
            double val = std::abs(std::round(cap(i, j) * n_nodes));
            g.add_edge(i, j, val, val);
        }
    }

    double flow = g.max_flow();
    return 0;
}*/

/*int main() {

    // std::vector<double> srcs_cap{0, 2};
    // std::vector<double> sink_cap{4, 7};
    // GraphCut g(2, srcs_cap.data(), sink_cap.data());
    // g.add_edge(0, 1, 5, 4);
    // const double flow = g.max_flow();
    std::vector<double> srcs_cap, sink_cap;
    std::vector<std::unordered_map<uint32_t, double>> edges;
    read_flow(srcs_cap, sink_cap, edges);
    GraphCut g(sink_cap.size(), srcs_cap.data(), sink_cap.data());
    for (uint32_t a = 0; a < edges.size(); a++) {
        for (const auto& pair : edges[a]) {
            g.add_edge(a, pair.first, pair.second, pair.second);
        }
    }
    const double flow = g.max_flow();
    return 0;
}*/

int main() {
    exactinit();
    // const uint32_t n_points = 10000;
    // auto points = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>::Random(n_points, 3).eval();
    // const auto tets = TetMesh::tetrahedralize(points.data(), n_points, 1e-6);
    const uint32_t n_points = 8;
    std::vector<double> points{
        0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1,
    };
    std::vector<uint32_t> indices{0, 3, 2, 0, 4, 3, 3, 4, 7, 3, 7, 6, 3, 6, 2, 6, 5,
                                  2, 2, 5, 1, 1, 4, 0, 1, 5, 4, 4, 5, 6, 4, 6, 7};
    make_polyhedral_mesh_from_triangles(points.data(), n_points, indices.data(), 11);
    return 0;
}
