#include "nlohmann/json.hpp"
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

uint32_t make_polyhedral_mesh(
    const double* points, const uint32_t n_points, const uint32_t* edge_data, const double* axis_data,
    const uint32_t* seperator, const uint32_t n_polygons, double** out_points, uint32_t** out_loops,
    uint32_t** out_loop_separators, uint32_t** out_polys, uint32_t** out_poly_separators, double** out_axis_data
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

static void writePolygon(
    const std::string& name, const double* data, const uint32_t n_points, const uint32_t* indices,
    const uint32_t* seperator, const uint32_t n_polygon
) {
    std::ofstream file(name);
    file.precision(std::numeric_limits<double>::max_digits10);
    for (uint32_t i = 0; i < n_points; i++) {
        file << "v " << data[3 * i] << " " << data[3 * i + 1] << " " << data[3 * i + 2] << "\n";
    }
    for (uint32_t i = 0; i < n_polygon; i++) {
        file << "f";
        for (uint32_t j = seperator[i]; j < seperator[i + 1]; j++) {
            file << " " << indices[j] + 1;
        }
        file << "\n";
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

template <typename scalar, typename Index>
inline bool
readOBJ(const std::string& file_name, std::vector<std::vector<scalar>>& V, std::vector<std::vector<Index>>& F) {
    FILE* obj_file = fopen(file_name.c_str(), "r");
    if (obj_file == NULL) {
        std::cout << "maybe error in filename\n";
        return false;
    }
    V.clear();
    F.clear();

    std::string v("v");
    std::string vn("vn");
    std::string vt("vt");
    std::string f("f");

#define LINE_MAX_ 2048
    char line[LINE_MAX_];
    std::vector<scalar> tv(3);
    int count = 0;
    while (fgets(line, LINE_MAX_, obj_file)) {
        char type[LINE_MAX_];
        if (sscanf(line, "%s", type) == 1) {
            char* c = &line[strlen(type)];
            if (type == v) {
                count = sscanf(c, "%lf %lf %lf\n", &tv[0], &tv[1], &tv[2]);
                if (count < 3) {
                    std::cout << "ERROR:vertex should have at least 3 coordinates\n";
                    fclose(obj_file);
                    return false;
                }
                V.push_back(tv);
            } else if (type == vn) {
                count = sscanf(c, "%lf %lf %lf\n", &tv[0], &tv[1], &tv[2]);
                if (count != 3) {
                    std::cout << "ERROR: normal should have 3 coordinates\n";
                    fclose(obj_file);
                    return false;
                }
            } else if (type == vt) {
                count = sscanf(c, "%lf %lf %lf\n", &tv[0], &tv[1], &tv[2]);
                if (count == 3) {
                } else if (count == 2) {
                } else {
                    std::cout << "ERROR:texture should have 3 or 2 coordinates\n";
                    fclose(obj_file);
                    return false;
                }
            } else if (type == f) {
                const auto& shift = [&V](const int i) -> int { return i < 0 ? i + V.size() : i - 1; };
                std::vector<Index> f;
                // Read each "word" after type
                char word[LINE_MAX_];
                int offset;
                while (sscanf(c, "%s%n", word, &offset) == 1) {
                    // adjust offset
                    c += offset;
                    // Process word
                    long int i, it, in;
                    if (sscanf(word, "%ld/%ld/%ld", &i, &it, &in) == 3) {
                        f.push_back(shift(i));

                    } else if (sscanf(word, "%ld/%ld", &i, &it) == 2) {
                        f.push_back(shift(i));
                    } else if (sscanf(word, "%ld//%ld", &i, &in) == 2) {
                        f.push_back(shift(i));
                    } else if (sscanf(word, "%ld", &i) == 1) {
                        f.push_back(shift(i));
                    } else {
                        std::cout << "Error: readOBJ() face on line has invalid element format\n";
                        fclose(obj_file);
                        return false;
                    }
                }
                if (f.size() > 0) {
                    // No matter what add each type to lists so that lists are the
                    // correct lengths
                    F.push_back(f);
                } else {
                    std::cout << "Error: readOBJ() face on line has invalid format\n";
                    fclose(obj_file);
                    return false;
                }

            } else if (strlen(type) >= 1 && (type[0] == '#' || type[0] == 'g' || type[0] == 's' || strcmp("usemtl", type) == 0 || strcmp("mtllib", type) == 0)) {
                // ignore comments or other shit
            } else {
                std::cout << "Warning: readOBJ() ignored non-comment line \n";
            }
        } else {
            // ignore empty line
        }
    }
    fclose(obj_file);

    return true;
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

static void read_input(
    const std::string& name, std::vector<double>& point_data, std::vector<uint32_t>& indices,
    std::vector<uint32_t>& separators, std::vector<double>& axes
) {
    using json = nlohmann::json;
    std::ifstream f(name);
    const json data = json::parse(f);
    data["pointsData"].get_to<std::vector<double>>(point_data);
    data["separators"].get_to<std::vector<uint32_t>>(separators);
    data["edgesData"].get_to<std::vector<uint32_t>>(indices);
    data["axesData"].get_to<std::vector<double>>(axes);
}

int main() {
    exactinit();
    std::vector<double> point_data;
    std::vector<uint32_t> indices;
    std::vector<uint32_t> seperator;
    std::vector<double> axes;
    read_input("input.json", point_data, indices, seperator, axes);
    double* out_points;
    double* out_axes;
    uint32_t* out_loops;
    uint32_t* out_loop_separators;
    uint32_t* out_polys;
    uint32_t* out_poly_separators;
    const auto n_polygons = make_polyhedral_mesh(
        point_data.data(), static_cast<uint32_t>(point_data.size() / 3), indices.data(), axes.data(), seperator.data(),
        static_cast<uint32_t>(seperator.size() - 1), &out_points, &out_loops, &out_loop_separators, &out_polys, &out_poly_separators, &out_axes
    );
    const auto out_n_pts = *std::max_element(out_loops, out_loops + out_loop_separators[n_polygons]) + 1;
    writePolygon("123.obj", out_points, out_n_pts, out_loops, out_loop_separators, n_polygons);
    delete[] out_points;
    delete[] out_axes;
    delete[] out_loops;
    delete[] out_loop_separators;
    delete[] out_polys;
    delete[] out_poly_separators;
}
