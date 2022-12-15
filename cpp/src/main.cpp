#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <fstream>
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

using SPMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using T = Eigen::Triplet<double>;


class GraphCut {
    static constexpr uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    struct Arc {
        uint32_t head;          // head
        uint32_t next{INVALID}; // next
        uint32_t sister;        // sister
        double r_cap {0};       // residual capacity
        Arc() {}
        Arc(const uint32_t h, const uint32_t n, const uint32_t s, const double c): head(h), next(n), sister(s), r_cap(c) {}
    };
        
    const uint32_t n_nodes;
    const uint32_t TERMINAL;
    const uint32_t ORPHAN;
    uint32_t time {0};
    uint32_t queue_first {INVALID};
    uint32_t queue_last {INVALID};
    std::vector<double> tr_cap;      // residual capacities of nodes
    double flow{0.0};
    std::vector<Arc> arcs;           // all arcs
    std::vector<uint32_t> dist;      // distance to the terminal
    std::vector<uint32_t> first_arc; // first arcs of nodes
    std::vector<uint32_t> next;      // next active node
    std::vector<uint32_t> parent;    // parent arcs of nodes
    std::vector<uint32_t> ts;        // time stamp showing when dist was computed
    std::vector<bool> is_sink;       // whether the node is in the source

  public:
    GraphCut(const uint32_t n, double* source_cap, double* sink_cap) : n_nodes(n), TERMINAL(n * 2), ORPHAN(TERMINAL + 1) {
        dist.resize(n, 1);
        is_sink.resize(n);
        next.resize(n, INVALID);
        parent.resize(n, TERMINAL);
        ts.resize(n, 0);
        for (uint32_t i = 0; i < n; i++) {
            tr_cap[i] = source_cap[i] - sink_cap[i];
            flow += source_cap[i] < sink_cap[i] ? source_cap[i] : sink_cap[i];
            is_sink[i] = tr_cap[i] < 0.0;
            if (tr_cap[i] == 0.0) {
                parent[i] = INVALID;
            } else {
                set_active(i);
            }
        }
        first_arc.resize(n, INVALID);
        
    }
        
    void set_active(const uint32_t i) {
        if (queue_last == INVALID) {
            queue_first = i;
        } else {
            next[queue_last] = i;
        }
        queue_last = i;
        next[i] = i;
    }
    
    uint32_t next_active() {
        while (true) {
            uint32_t i = queue_first;
            if (next[i] == i) {
                queue_first = queue_last = INVALID;
            } else {
                queue_first = next[i];
            }
            next[i] = INVALID;
            if (parent[i] != INVALID) {
                return i;
            }
        }
    }
        
    void add_edge(const uint32_t i, const uint32_t j, const double cap, const double rev_cap) {
        const uint32_t arc_id = static_cast<uint32_t>(arcs.size());
        const uint32_t sisiter_arc_id = arc_id + 1;
        arcs.emplace_back(j, first_arc[i], sisiter_arc_id, cap); // arc
        arcs.emplace_back(i, first_arc[j], arc_id, rev_cap);     // sisiter arc
    }
    double max_flow() {
        uint32_t current_node = INVALID;
        while (true) {
            uint32_t i = current_node;
            if (i != INVALID) {
                next[i] = INVALID;
                if (parent[i] == INVALID) {
                    i = INVALID;
                }
            }
            if (i == INVALID) {
                if ((i = next_active()) == INVALID) {
                    break;
                }
            }
            uint32_t aid = INVALID;
            if (!is_sink[i]) {
                aid = first_arc[i];
                while (aid != INVALID) {
                    const Arc& a = arcs[aid];
                    if (a.r_cap != 0.0) {
                        const uint32_t j = a.head;
                        if (parent[j] == INVALID) {
                            is_sink[j] = false;
                            parent[j] = a.sister;
                            ts[j] = ts[i];
                            set_active(j);
                        } else if (is_sink[j]) {
                            break;
                        } else if (ts[j] <= ts[i] && dist[j] > dist[i]) {
                            parent[j] = a.sister;
                            ts[j] = ts[i];
                            dist[j] = dist[i] + 1;
                        }
                    }
                    aid = a.next;
                }
            } else {
                aid = first_arc[i];
                while (aid != INVALID) {
                    const Arc& a = arcs[aid];
                    if (a.r_cap != 0.0) {
                        const uint32_t j = a.head;
                        if (parent[j] == INVALID) {
                            is_sink[j] = true;
                            parent[j] = a.sister;
                            ts[j] = ts[i];
                        } else if (!is_sink[j]) {
                            break;
                        } else if (ts[j] <= ts[i] && dist[j] > dist[i]) {
                            parent[j] = a.sister;
                            ts[j] = ts[i];
                            dist[j] = dist[i + 1];
                                        
                        }
                    }
                    aid = a.next;
                }
            }
            time += 1;
            if (aid != INVALID) {
            
            }
        }
        return flow;
    }
};

int main() {
    std::vector<T> triples{{1, 0, 22.0}, {2, 0, 7.0}, {0, 1, 3.0},  {2, 1, 5.0},
                           {4, 2, 14.0}, {2, 3, 2.0}, {1, 4, 17.0}, {4, 4, 8.0}};
    const uint32_t n = 5;
    SPMat m(n, n);
    m.setFromTriplets(triples.begin(), triples.end());
    m.coeffRef(0, 1) = 0.0;
    m.prune(0.0);
    for (int k = 0; k < m.outerSize(); ++k) {
        std::cout << "k: " << k << "\n";
        for (SPMat::InnerIterator it(m, k); it; ++it) {
            it.value();
            it.row();   // row index
            it.col();   // col index (here it is equal to k)
            it.index(); // inner index, here it is equal to it.row()
            std::cout << "(" << it.row() << ", " << it.col() << "): " << it.value() << "\n";
        }
    }
    return 0;
}

/*int main() {
    exactinit();
    // const uint32_t n_points = 10000;
    // auto points = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>::Random(n_points, 3).eval();
    // const auto tets = TetMesh::tetrahedralize(points.data(), n_points, 1e-6);
    const uint32_t n_points = 8;
    std::vector<double> points {
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        0, 0, 1,
        1, 0, 1,
        1, 1, 1,
        0, 1, 1,
    };
    std::vector<uint32_t> indices{0, 3, 2, 0, 4, 3, 3, 4, 7, 3, 7, 6, 3, 6, 2, 6, 5, 2, 2, 5, 1, 1, 4, 0, 1, 5, 4, 4, 5,
6, 4, 6, 7}; make_polyhedral_mesh_from_triangles(points.data(), n_points, indices.data(), 11); return 0;
}*/
