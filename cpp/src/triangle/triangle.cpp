#include <triangle/triangle.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

extern "C" {
double orient2d(double*, double*, double*);
double incircle(double* pa, double* pb, double* pc, double* pd);
}

constexpr uint32_t INVALID{std::numeric_limits<uint32_t>::max()};
struct HEdge {
    uint32_t tri{INVALID};
    uint32_t ori{0};

    inline friend bool operator==(const HEdge& a, const HEdge& b) { return a.tri == b.tri && a.ori == b.ori; }
    inline friend bool operator!=(const HEdge& a, const HEdge& b) { return !(a == b); }
};

struct Triangle {
    std::array<uint32_t, 6> data{{INVALID, INVALID, INVALID, INVALID, INVALID, INVALID}};
    std::array<HEdge, 3> nei;
};

struct Mesh {
    const double* points;
    std::vector<Triangle> triangles;
};

static void alternate_axes(const double* points, uint32_t* indices, const uint32_t len, uint32_t axis) {
    uint32_t divider = len >> 1;
    if (len <= 3) {
        axis = 0;
    }
    std::nth_element(indices, indices + divider, indices + len, [&points, axis](const uint32_t i, const uint32_t j) {
        const double* pi = &points[i << 1];
        const double* pj = &points[j << 1];
        return pi[axis] < pj[axis] || (pi[axis] == pj[axis] && pi[1 - axis] < pj[1 - axis]);
    });
    if (len - divider >= 2) {
        if (divider >= 2) {
            alternate_axes(points, &indices[0], divider, 1 - axis);
        }
        alternate_axes(points, &indices[divider], len - divider, 1 - axis);
    }
}

inline void make_triangle(std::vector<Triangle>& triangles, HEdge& he) {
    he.tri = static_cast<uint32_t>(triangles.size());
    triangles.emplace_back();
    he.ori = 0;
}

inline void set_org(std::vector<Triangle>& triangles, const HEdge& he, const uint32_t vid) {
    triangles[he.tri].data[(he.ori + 1) % 3] = vid;
}

inline void set_dest(std::vector<Triangle>& triangles, const HEdge& he, const uint32_t vid) {
    triangles[he.tri].data[(he.ori + 2) % 3] = vid;
}

inline void set_apex(std::vector<Triangle>& triangles, const HEdge& he, const uint32_t vid) {
    triangles[he.tri].data[he.ori] = vid;
}

inline void bond(std::vector<Triangle>& triangles, HEdge& he1, HEdge& he2) {
    triangles[he1.tri].nei[he1.ori] = he2;
    triangles[he2.tri].nei[he2.ori] = he1;
}

inline uint32_t org(const Mesh* m, const HEdge& he) { return m->triangles[he.tri].data[(he.ori + 1) % 3]; }

inline uint32_t dest(const Mesh* m, const HEdge& he) { return m->triangles[he.tri].data[(he.ori + 2) % 3]; }

inline uint32_t apex(const Mesh* m, const HEdge& he) { return m->triangles[he.tri].data[he.ori]; }

inline const double* point(const double* points, const uint32_t idx) { return &points[idx << 1]; }

inline void prev_self(HEdge& he) { he.ori = (he.ori + 2) % 3; }

inline void next_self(HEdge& he) { he.ori = (he.ori + 1) % 3; }

inline void sym_self(const std::vector<Triangle>& triangles, HEdge& he) {
    const auto& nei = triangles[he.tri].nei[he.ori];
    he.tri = nei.tri;
    he.ori = nei.ori;
}

inline void oprev_self(const std::vector<Triangle>& triangles, HEdge& he) {
    sym_self(triangles, he);
    next_self(he);
}

inline void onext_self(const std::vector<Triangle>& triangles, HEdge& he) {
    prev_self(he);
    sym_self(triangles, he);
}


inline void prev(const HEdge& src, HEdge& dest) {
    dest.tri = src.tri;
    dest.ori = (src.ori + 2) % 3;
}

inline void next(const HEdge& src, HEdge& dest) {
    dest.tri = src.tri;
    dest.ori = (src.ori + 1) % 3;
}

inline void sym(const std::vector<Triangle>& triangles, const HEdge& src, HEdge& dest) {
    const auto& nei = triangles[src.tri].nei[src.ori];
    dest.tri = nei.tri;
    dest.ori = nei.ori;
}

inline void oprev(const std::vector<Triangle>& triangles, const HEdge& src, HEdge& dest) {
    sym(triangles, src, dest);
    next_self(dest);
}

inline void onext(const std::vector<Triangle>& triangles, const HEdge& src, HEdge& dest) {
    prev(src, dest);
    sym_self(triangles, dest);
}

inline void copy(const HEdge& src, HEdge& dest) {
    dest.tri = src.tri;
    dest.ori = src.ori;
}

inline double counterclockwise(const double* points, const uint32_t i, const uint32_t j, const uint32_t k) {
    return orient2d(
        const_cast<double*>(&points[i << 1]), const_cast<double*>(&points[j << 1]), const_cast<double*>(&points[k << 1])
    );
}

inline double incircle(const double* points, const uint32_t a, const uint32_t b, const uint32_t c, const uint32_t d) {
    return incircle(
        const_cast<double*>(&points[a << 1]), const_cast<double*>(&points[b << 1]),
        const_cast<double*>(&points[c << 1]), const_cast<double*>(&points[d << 1])
    );
}

inline void
merge_hulls(Mesh* m, const uint32_t axis, HEdge& far_left, HEdge& inner_left, HEdge& inner_right, HEdge& far_right) {
    uint32_t inner_left_dest = dest(m, inner_left);
    uint32_t inner_left_apex = apex(m, inner_left);
    uint32_t inner_right_org = org(m, inner_right);
    uint32_t inner_right_apex = apex(m, inner_right);
    uint32_t far_left_pt, far_left_apex, far_right_pt, far_right_apex;
    // Special treatment for horizontal cuts
    if (axis == 1) {
        far_left_pt = org(m, far_left);
        far_left_apex = apex(m, far_left);
        far_right_pt = dest(m, far_right);
        far_right_apex = apex(m, far_right);
        // The pointers to the extremal vertices are shifted to point to the
        // topmost and bottommost vertex of each hull, rather than the
        // leftmost and rightmost vertices.
        while (point(m->points, far_left_apex)[1] < point(m->points, far_left_pt)[1]) {
            next_self(far_left);
            sym_self(m->triangles, far_left);
            far_left_pt = far_left_apex;
            far_left_apex = apex(m, far_left);
        }

        HEdge check_edge;
        sym(m->triangles, inner_left, check_edge);
        uint32_t check_vertex = apex(m, check_edge);
        while (point(m->points, check_vertex)[1] > point(m->points, inner_left_dest)[1]) {
            next(check_edge, inner_left);
            inner_left_apex = inner_left_dest;
            inner_left_dest = check_vertex;
            sym(m->triangles, inner_left, check_edge);
            check_vertex = apex(m, check_edge);
        }

        while (point(m->points, inner_right_apex)[1] < point(m->points, inner_right_org)[1]) {
            next_self(inner_right);
            sym_self(m->triangles, inner_right);
            inner_right_org = inner_right_apex;
            inner_right_apex = apex(m, inner_right);
        }

        sym(m->triangles, far_right, check_edge);
        check_vertex = apex(m, check_edge);
        while (point(m->points, check_vertex)[1] > point(m->points, far_right_pt)[1]) {
            next(check_edge, far_right);
            far_right_apex = far_right_pt;
            far_right_pt = check_vertex;
            sym(m->triangles, far_right, check_edge);
            check_vertex = apex(m, check_edge);
        }
    }

    // Find a line tangent to and below both hulls
    bool change_made = true;
    while (change_made) {
        change_made = false;
        // Make innerleftdest the "bottommost" vertex of the left hull
        if (counterclockwise(m->points, inner_left_dest, inner_left_apex, inner_right_org) > 0.0) {
            prev_self(inner_left);
            sym_self(m->triangles, inner_left);
            inner_left_dest = inner_left_apex;
            inner_left_apex = apex(m, inner_left);
            change_made = true;
        }
        // Make innerrightorg the "bottommost" vertex of the right hull
        if (counterclockwise(m->points, inner_right_apex, inner_right_org, inner_left_dest) > 0.0) {
            next_self(inner_right);
            sym_self(m->triangles, inner_right);
            inner_right_org = inner_right_apex;
            inner_right_apex = apex(m, inner_right);
            change_made = true;
        }
    }

    // Find the two candidates to be the next "gear tooth"
    HEdge left_cand, right_cand, base_edge;
    sym(m->triangles, inner_left, left_cand);
    sym(m->triangles, inner_right, right_cand);
    // Create the bottom new bounding triangle
    make_triangle(m->triangles, base_edge);
    // Connect it to the bounding boxes of the left and right triangulations.
    bond(m->triangles, base_edge, inner_left);
    next_self(base_edge);
    bond(m->triangles, base_edge, inner_right);
    next_self(base_edge);
    set_org(m->triangles, base_edge, inner_right_org);
    set_dest(m->triangles, base_edge, inner_left_dest);


    // Fix the extreme triangles if necessary
    far_left_pt = org(m, far_left);
    if (inner_left_dest == far_left_pt) {
        next(base_edge, far_left);
    }
    far_right_pt = dest(m, far_right);
    if (inner_right_org == far_right_pt) {
        prev(base_edge, far_right);
    }

    // The vertices of the current knitting edge
    uint32_t lower_left = inner_left_dest;
    uint32_t lower_right = inner_right_org;
    // The candidate vertices for knitting
    uint32_t upper_left = apex(m, left_cand);
    uint32_t upper_right = apex(m, right_cand);
    // Walk up the gap between the two triangulations, knitting them together
    while (true) {
        const auto left_finished = counterclockwise(m->points, upper_left, lower_left, lower_right) <= 0.0;
        const auto right_finished = counterclockwise(m->points, upper_right, lower_left, lower_right) <= 0.0;
        HEdge check_edge, next_edge;
        if (left_finished && right_finished) {
            // Create the top new bounding triangle
            make_triangle(m->triangles, next_edge);
            set_org(m->triangles, next_edge, lower_left);
            set_dest(m->triangles, next_edge, lower_right);
            // Apex is intentionally left INVALID
            // Connect it to the bounding boxes of the two triangulations
            bond(m->triangles, next_edge, base_edge);
            next_self(next_edge);
            bond(m->triangles, next_edge, right_cand);
            next_self(next_edge);
            bond(m->triangles, next_edge, left_cand);

            // Special treatment for horizontal cuts
            if (axis == 1) {
                far_left_pt = org(m, far_left);
                far_left_apex = apex(m, far_left);
                far_right_pt = dest(m, far_right);
                far_right_apex = apex(m, far_right);
                sym(m->triangles, far_left, check_edge);
                uint32_t check_vertex = apex(m, check_edge);
                // The pointers to the extremal vertices are restored to the
                // leftmost and rightmost vertices (rather than topmost and
                // bottommost)
                while (point(m->points, check_vertex)[0] < point(m->points, far_left_pt)[0]) {
                    prev(check_edge, far_left);
                    far_left_apex = far_left_pt;
                    far_left_pt = check_vertex;
                    sym(m->triangles, far_left, check_edge);
                    check_vertex = apex(m, check_edge);
                }
                while (point(m->points, far_right_apex)[0] > point(m->points, far_right_pt)[0]) {
                    prev_self(far_right);
                    sym_self(m->triangles, far_right);
                    far_right_pt = far_right_apex;
                    far_right_apex = apex(m, far_right);
                }
            }
            return;
        }

        // Consider eliminating edges from the left triangulation
        if (!left_finished) {
            // What vertex would be exposed if an edge were deleted
            prev(left_cand, next_edge);
            sym_self(m->triangles, next_edge);
            uint32_t next_apex = apex(m, next_edge);
            // If nextapex is INVALID, then no vertex would be exposed; the
            // triangulation would have been eaten right through.
            if (next_apex != INVALID) {
                // Check whether the edge is Delaunay
                auto bad_edge = incircle(m->points, lower_left, lower_right, upper_left, next_apex) > 0.0;
                HEdge top_casing, side_casing, outer_casing;
                while (bad_edge) {
                    // Eliminate the edge with an edge flip.  As a result, the
                    // left triangulation will have one more boundary triangle.
                    next_self(next_edge);
                    sym(m->triangles, next_edge, top_casing);
                    next_self(next_edge);
                    sym(m->triangles, next_edge, side_casing);
                    bond(m->triangles, next_edge, top_casing);
                    bond(m->triangles, left_cand, side_casing);
                    next_self(left_cand);
                    sym(m->triangles, left_cand, outer_casing);
                    prev_self(next_edge);
                    bond(m->triangles, next_edge, outer_casing);
                    // Correct the vertices to reflect the edge flip
                    set_org(m->triangles, left_cand, lower_left);
                    set_dest(m->triangles, left_cand, INVALID);
                    set_apex(m->triangles, left_cand, next_apex);
                    set_org(m->triangles, next_edge, INVALID);
                    set_dest(m->triangles, next_edge, upper_left);
                    set_apex(m->triangles, next_edge, next_apex);
                    // Consider the newly exposed vertex
                    upper_left = next_apex;
                    // What vertex would be exposed if another edge were deleted?
                    copy(side_casing, next_edge);
                    next_apex = apex(m, next_edge);
                    if (next_apex != INVALID) {
                        // Check whether the edge is Delaunay
                        bad_edge = incircle(m->points, lower_left, lower_right, upper_left, next_apex) > 0.0;
                    } else {
                        // Avoid eating right through the triangulation
                        bad_edge = false;
                    }
                }
            }
        }

        // Consider eliminating edges from the right triangulation
        if (!right_finished) {
            next(right_cand, next_edge);
            sym_self(m->triangles, next_edge);
            uint32_t next_apex = apex(m, next_edge);
            if (next_apex != INVALID) {
                // Check whether the edge is Delaunay
                bool bad_edge = incircle(m->points, lower_left, lower_right, upper_right, next_apex) > 0.0;
                HEdge top_casing, side_casing, outer_casing;
                while (bad_edge) {
                    prev_self(next_edge);
                    sym(m->triangles, next_edge, top_casing);
                    prev_self(next_edge);
                    sym(m->triangles, next_edge, side_casing);
                    bond(m->triangles, next_edge, top_casing);
                    bond(m->triangles, right_cand, side_casing);
                    prev_self(right_cand);
                    sym(m->triangles, right_cand, outer_casing);
                    next_self(next_edge);
                    bond(m->triangles, next_edge, outer_casing);

                    set_org(m->triangles, right_cand, INVALID);
                    set_dest(m->triangles, right_cand, lower_right);
                    set_apex(m->triangles, right_cand, next_apex);
                    set_org(m->triangles, next_edge, upper_right);
                    set_dest(m->triangles, next_edge, INVALID);
                    set_apex(m->triangles, next_edge, next_apex);

                    upper_right = next_apex;

                    copy(side_casing, next_edge);
                    next_apex = apex(m, next_edge);
                    if (next_apex != INVALID) {
                        bad_edge = incircle(m->points, lower_left, lower_right, upper_right, next_apex) > 0.0;
                    } else {
                        bad_edge = false;
                    }
                }
            }
        }


        if (left_finished ||
            (!right_finished && (incircle(m->points, upper_left, lower_left, lower_right, upper_right) > 0.0))) {
            // Knit the triangulations, adding an edge from `lowerleft'
            // to `upperright'
            bond(m->triangles, base_edge, right_cand);
            prev(right_cand, base_edge);
            set_dest(m->triangles, base_edge, lower_left);
            lower_right = upper_right;
            sym(m->triangles, base_edge, right_cand);
            upper_right = apex(m, right_cand);
        } else {
            // Knit the triangulations, adding an edge from `upperleft'
            // to `lowerright'
            bond(m->triangles, base_edge, left_cand);
            next(left_cand, base_edge);
            set_org(m->triangles, base_edge, lower_right);
            lower_left = upper_left;
            sym(m->triangles, base_edge, left_cand);
            upper_left = apex(m, left_cand);
        }
    }
}

static void div_conq_recurse(
    Mesh* m, const std::vector<uint32_t>& sorted_pt_inds, const uint32_t axis, const uint32_t start, const uint32_t end,
    HEdge& far_left, HEdge& far_right
) {
    const auto len = end - start;
    if (len == 2) {
        make_triangle(m->triangles, far_left);
        set_org(m->triangles, far_left, sorted_pt_inds[start]);
        set_dest(m->triangles, far_left, sorted_pt_inds[start + 1]);

        make_triangle(m->triangles, far_right);
        set_org(m->triangles, far_right, sorted_pt_inds[start + 1]);
        set_dest(m->triangles, far_right, sorted_pt_inds[start]);
        bond(m->triangles, far_left, far_right);

        prev_self(far_left);
        next_self(far_right);
        bond(m->triangles, far_left, far_right);

        prev_self(far_left);
        next_self(far_right);
        bond(m->triangles, far_left, far_right);

        // ensura that the origin of `farleft` is `start`
        prev(far_right, far_left);
    } else if (len == 3) {
        HEdge midtri, tri1, tri2, tri3;
        make_triangle(m->triangles, midtri);
        make_triangle(m->triangles, tri1);
        make_triangle(m->triangles, tri2);
        make_triangle(m->triangles, tri3);
        const auto area =
            counterclockwise(m->points, sorted_pt_inds[start], sorted_pt_inds[start + 1], sorted_pt_inds[start + 2]);
        if (area == 0.0) {
            // Three collinear vertices; the triangulation is two edges
            set_org(m->triangles, midtri, sorted_pt_inds[start]);
            set_dest(m->triangles, midtri, sorted_pt_inds[start + 1]);
            set_org(m->triangles, tri1, sorted_pt_inds[start + 1]);
            set_dest(m->triangles, tri1, sorted_pt_inds[start]);
            set_org(m->triangles, tri2, sorted_pt_inds[start + 2]);
            set_dest(m->triangles, tri2, sorted_pt_inds[start + 1]);
            set_org(m->triangles, tri3, sorted_pt_inds[start + 1]);
            set_dest(m->triangles, tri3, sorted_pt_inds[start + 2]);
            // All apices are intentionally left INVALID.
            bond(m->triangles, midtri, tri1);
            bond(m->triangles, tri2, tri3);
            next_self(midtri);
            prev_self(tri1);
            next_self(tri2);
            prev_self(tri3);
            bond(m->triangles, midtri, tri3);
            bond(m->triangles, tri1, tri2);
            next_self(midtri);
            prev_self(tri1);
            next_self(tri2);
            prev_self(tri3);
            bond(m->triangles, midtri, tri1);
            bond(m->triangles, tri2, tri3);

            copy(tri1, far_left);
            copy(tri2, far_right);
        } else {

            // The three vertices are not collinear; the triangulation is one
            // triangle, namely `midtri'.
            set_org(m->triangles, midtri, sorted_pt_inds[start]);
            set_dest(m->triangles, tri1, sorted_pt_inds[start]);
            set_org(m->triangles, tri3, sorted_pt_inds[start]);
            if (area > 0.0) {
                // The vertices are in counterclockwise order
                set_dest(m->triangles, midtri, sorted_pt_inds[start + 1]);
                set_org(m->triangles, tri1, sorted_pt_inds[start + 1]);
                set_dest(m->triangles, tri2, sorted_pt_inds[start + 1]);
                set_apex(m->triangles, midtri, sorted_pt_inds[start + 2]);
                set_org(m->triangles, tri2, sorted_pt_inds[start + 2]);
                set_dest(m->triangles, tri3, sorted_pt_inds[start + 2]);
            } else {
                // The vertices are in clockwise order
                set_dest(m->triangles, midtri, sorted_pt_inds[start + 2]);
                set_org(m->triangles, tri1, sorted_pt_inds[start + 2]);
                set_dest(m->triangles, tri2, sorted_pt_inds[start + 2]);
                set_apex(m->triangles, midtri, sorted_pt_inds[start + 1]);
                set_org(m->triangles, tri2, sorted_pt_inds[start + 1]);
                set_dest(m->triangles, tri3, sorted_pt_inds[start + 1]);
            }
            // The topology does not depend on how the vertices are ordered
            bond(m->triangles, midtri, tri1);
            next_self(midtri);
            bond(m->triangles, midtri, tri2);
            next_self(midtri);
            bond(m->triangles, midtri, tri3);
            prev_self(tri1);
            next_self(tri2);
            bond(m->triangles, tri1, tri2);
            prev_self(tri1);
            prev_self(tri3);
            bond(m->triangles, tri1, tri3);
            next_self(tri2);
            prev_self(tri3);
            bond(m->triangles, tri2, tri3);
            // Ensure that the origin of `farleft' is start
            copy(tri1, far_left);
            // Ensure that the destination of `farright' is `start + 2`
            if (area > 0.0) {
                copy(tri2, far_right);
            } else {
                next(far_left, far_right);
            }
        }
    } else {
        const auto divider = len >> 1;
        HEdge innerleft, innerright;

        div_conq_recurse(m, sorted_pt_inds, 1 - axis, start, start + divider, far_left, innerleft);
        div_conq_recurse(m, sorted_pt_inds, 1 - axis, start + divider, end, innerright, far_right);
        merge_hulls(m, axis, far_left, innerleft, innerright, far_right);
    }
}

inline uint32_t mark_ghost(std::vector<Triangle>& triangles, const HEdge& start, std::vector<bool>& ghost) {
    HEdge dissolve_edge;
    copy(start, dissolve_edge);
    uint32_t count = 0;
    do {
        ghost[dissolve_edge.tri] = true;
        count += 1;
        next_self(dissolve_edge);
        sym_self(triangles, dissolve_edge);
    } while (dissolve_edge != start);
    return count;
}

inline std::vector<HEdge>
make_vertex_map(const std::vector<Triangle>& triangles, const std::vector<bool>& ghost, const uint32_t n_points) {
    std::vector<HEdge> vertex_map(n_points);
    for (uint32_t idx = 0; idx < triangles.size(); idx++) {
        if (ghost[idx]) {
            continue;
        }
        const auto& tri = triangles[idx];
        for (uint32_t i = 0; i < 3; i++) {
            const auto org = tri.data[i];
            if (vertex_map[org].tri != INVALID) {
                continue;
            }
            const auto he = tri.nei[(i + 1) % 3];
            if (!ghost[he.tri]) {
                vertex_map[org] = he;
            } else {
                vertex_map[org].tri = idx;
                vertex_map[org].ori = (i + 2) % 3;
            }
        }
    }
    return vertex_map;
}

enum class Direction { WITHIN, LEFTCOLLINEAR, RIGHTCOLLINEAR };

Direction find_direction(Mesh* m, HEdge& search_tri, uint32_t search_point) {

    uint32_t start_vertex = org(m, search_tri);
    uint32_t right_vertex = dest(m, search_tri);
    uint32_t left_vertex = apex(m, search_tri);

    double left_ccw = counterclockwise(m->points, search_point, start_vertex, left_vertex);
    bool left_flag = left_ccw > 0.0;

    double right_ccw = counterclockwise(m->points, start_vertex, search_point, right_vertex);
    bool right_flag = right_ccw > 0.0;
    HEdge check_tri;
    if (left_flag && right_flag) {
        onext(m->triangles, search_tri, check_tri);
        if (check_tri.tri == INVALID) {
            left_flag = false;
        } else {
            right_flag = false;
        }
    }

    uint32_t iter = 0;
    // Valid worst case: vertex is incident on every face.
    const uint32_t max_iter = static_cast<uint32_t>(m->triangles.size() << 1) + 100;
    while (left_flag) {
        // Turn left until satisfied
        onext_self(m->triangles, search_tri);
        if (search_tri.tri == INVALID) {
            throw "unreachable";
        }
        left_vertex = apex(m, search_tri);
        right_ccw = left_ccw;
        left_ccw = counterclockwise(m->points, search_point, start_vertex, left_vertex);
        left_flag = left_ccw > 0.0;
        iter++;
        if (iter > max_iter) {
            throw "infinite loop";
        }
    }

    iter = 0;
    while (right_flag) {
        // Turn right until satisfied
        oprev_self(m->triangles, search_tri);
        if (search_tri.tri == INVALID) {
            throw "unreachable";
        }
        right_vertex = dest(m, search_tri);
        left_ccw = right_ccw;
        right_ccw = counterclockwise(m->points, start_vertex, search_point, right_vertex);
        right_flag = right_ccw > 0.0;
        iter++;
        if (iter > max_iter) {
            throw "infinite loop";
        }
    }
    if (left_ccw == 0.0) {
        return Direction::LEFTCOLLINEAR;
    } else if (right_ccw == 0.0) {
        return Direction::RIGHTCOLLINEAR;
    } else {
        return Direction::WITHIN;
    }
}

inline uint32_t twin(const uint32_t mark) { return ((mark & 1) == 0) ? (mark + 1) : (mark - 1); }

inline void set_mark(std::vector<Triangle>& triangles, const HEdge& he, const uint32_t mark) {
    triangles[he.tri].data[he.ori + 3] = mark;
    HEdge sym_edge;
    sym(triangles, he, sym_edge);
    if (sym_edge.tri != INVALID) {
        triangles[sym_edge.tri].data[sym_edge.ori + 3] = twin(mark);
    }
}

inline uint32_t mark(std::vector<Triangle>& triangles, const HEdge& he) { return triangles[he.tri].data[he.ori + 3]; }

static bool scout_segment(Mesh* m, HEdge& search_tri, const uint32_t endpoint2, const uint32_t mark) {

    const auto collinear = find_direction(m, search_tri, endpoint2);
    const uint32_t right_vertex = dest(m, search_tri);
    const uint32_t left_vertex = apex(m, search_tri);
    if (left_vertex == endpoint2) {
        // The segment is already an edge in the mesh.
        prev_self(search_tri);
        set_mark(m->triangles, search_tri, twin(mark));
        return true;
    } else if (right_vertex == endpoint2) {
        // The segment is already an edge in the mesh.
        set_mark(m->triangles, search_tri, mark);
        return true;
    } else if (collinear == Direction::LEFTCOLLINEAR) {
        prev_self(search_tri);
        return scout_segment(m, search_tri, endpoint2, twin(mark));
    } else if (collinear == Direction::RIGHTCOLLINEAR) {
        next_self(search_tri);
        return scout_segment(m, search_tri, endpoint2, mark);
    } else {
        return false;
    }
}

inline void flip(Mesh* m, HEdge& flip_edge) {
    uint32_t right_vertex = org(m, flip_edge);
    uint32_t left_vertex = dest(m, flip_edge);
    uint32_t bot_vertex = apex(m, flip_edge);
    HEdge top;
    sym(m->triangles, flip_edge, top);
    uint32_t far_vertex = apex(m, top);

    // Identify the casing of the quadrilateral.
    HEdge top_left;
    prev(top, top_left);
    HEdge topl_casing;
    sym(m->triangles, top_left, topl_casing);
    HEdge top_right;
    next(top, top_right);
    HEdge topr_casing;
    sym(m->triangles, top_right, topr_casing);
    HEdge bot_left;
    next(flip_edge, bot_left);
    HEdge botl_casing;
    sym(m->triangles, bot_left, botl_casing);
    HEdge bot_right;
    prev(flip_edge, bot_right);
    HEdge botr_casing;
    sym(m->triangles, bot_right, botr_casing);
    // Rotate the quadrilateral one-quarter turn counterclockwise
    bond(m->triangles, top_left, botl_casing);
    bond(m->triangles, bot_left, botr_casing);
    bond(m->triangles, bot_right, topr_casing);
    bond(m->triangles, top_right, topl_casing);

    set_org(m->triangles, flip_edge, far_vertex);
    set_dest(m->triangles, flip_edge, bot_vertex);
    set_apex(m->triangles, flip_edge, right_vertex);
    set_org(m->triangles, top, bot_vertex);
    set_dest(m->triangles, top, far_vertex);
    set_apex(m->triangles, top, left_vertex);
}

static void delaunay_fixup(Mesh* m, HEdge& fixup_tri, const bool left_side) {
    HEdge near_tri;
    next(fixup_tri, near_tri);
    HEdge far_tri;
    sym(m->triangles, near_tri, far_tri);
    // Check if the edge opposite the origin of fixuptri can be flipped
    if (far_tri.tri == INVALID) {
        return;
    }
    if (mark(m->triangles, near_tri) != INVALID) {
        return;
    }
    uint32_t near_vertex = apex(m, near_tri);
    uint32_t left_vertex = org(m, near_tri);
    uint32_t right_vertex = dest(m, near_tri);
    uint32_t far_vertex = apex(m, far_tri);
    // Check whether the previous polygon vertex is a reflex vertex
    if (left_side) {
        if (counterclockwise(m->points, near_vertex, left_vertex, far_vertex) <= 0.0) {
            // leftvertex is a reflex vertex too, nothing can
            // be done until a convex section is found.
            return;
        }
    } else {
        if (counterclockwise(m->points, far_vertex, right_vertex, near_vertex) <= 0.0) {
            return;
        }
    }

    if (counterclockwise(m->points, right_vertex, left_vertex, far_vertex) > 0.0) {
        if (incircle(m->points, left_vertex, far_vertex, right_vertex, near_vertex) <= 0.0) {
            return;
        }
    }

    flip(m, near_tri);
    prev_self(fixup_tri);
    delaunay_fixup(m, fixup_tri, left_side);
    delaunay_fixup(m, far_tri, left_side);
}

static void constrained_edge(Mesh* m, HEdge& start_tri, const uint32_t endpoint2, const uint32_t mark) {
    const auto endpoint1 = org(m, start_tri);
    HEdge fixup_tri;
    next(start_tri, fixup_tri);
    flip(m, fixup_tri);

    bool collision = 0;
    bool done = 0;
    do {
        uint32_t far_vertex = org(m, fixup_tri);
        if (far_vertex == endpoint2) {
            HEdge fixup_tri2;
            oprev(m->triangles, fixup_tri, fixup_tri2);
            // Enforce the Delaunay condition around endpoint2
            delaunay_fixup(m, fixup_tri, false);
            delaunay_fixup(m, fixup_tri2, true);
            done = true;
        } else {
            const double area = counterclockwise(m->points, endpoint1, endpoint2, far_vertex);
            if (area == 0.0) {
                collision = true;
                HEdge fixup_tri2;
                oprev(m->triangles, fixup_tri, fixup_tri2);

                delaunay_fixup(m, fixup_tri, false);
                delaunay_fixup(m, fixup_tri2, true);

                done = 1;
            } else {
                if (area > 0.0) {
                    HEdge fixuptri2;
                    oprev(m->triangles, fixup_tri, fixuptri2);
                    delaunay_fixup(m, fixuptri2, true);
                    prev_self(fixup_tri);
                } else {
                    delaunay_fixup(m, fixup_tri, false);
                    oprev_self(m->triangles, fixup_tri);
                }
                flip(m, fixup_tri);
            }
        }
    } while (!done);
    if (collision) {
        if (!scout_segment(m, fixup_tri, endpoint2, mark)) {
            constrained_edge(m, fixup_tri, endpoint2, mark);
        }
    } else {
        set_mark(m->triangles, fixup_tri, mark);
    }
}

inline void
insert_segment(Mesh* m, const std::vector<HEdge> vertex_map, uint32_t start, uint32_t end, const uint32_t mark) {
    HEdge searchtri1 = vertex_map[start];
    if (scout_segment(m, searchtri1, end, mark)) {
        return;
    }

    // The first endpoint may have changed if a collision with an intervening
    // vertex on the segment occurred.
    start = org(m, searchtri1);

    HEdge searchtri2 = vertex_map[end];
    if (scout_segment(m, searchtri2, start, twin(mark))) {
        return;
    }
    end = org(m, searchtri2);
    constrained_edge(m, searchtri1, end, mark);
}

static void form_skeleton(
    Mesh* m, const uint32_t n_points, const std::vector<bool>& ghost, const uint32_t* segments,
    const uint32_t n_segments
) {
    const auto vertex_map = make_vertex_map(m->triangles, ghost, n_points);
    for (uint32_t i = 0; i < n_segments; i++) {
        const auto* data = &segments[i << 1];
        if (data[0] != data[1]) {
            insert_segment(m, vertex_map, data[0], data[1], i << 1);
        }
    }
}

#ifdef EMSCRIPTEN
extern "C" {
#endif
uint32_t* triangulate(
    const double* points, uint32_t n_points, const uint32_t* segments, const uint32_t n_segments, uint32_t* n_triangles
) {
    std::vector<uint32_t> sorted_pt_inds(n_points);
    std::iota(sorted_pt_inds.begin(), sorted_pt_inds.end(), 0);
    std::sort(sorted_pt_inds.begin(), sorted_pt_inds.end(), [points](const uint32_t i, const uint32_t j) {
        const double* pi = &points[i << 1];
        const double* pj = &points[j << 1];
        return pi[0] < pj[0] || (pi[0] == pj[0] && pi[1] < pj[1]);
    });
    {
        /** todo: remove*/
        const auto it =
            std::unique(sorted_pt_inds.begin(), sorted_pt_inds.end(), [points](const uint32_t i, const uint32_t j) {
                const double* pi = &points[i << 1];
                const double* pj = &points[j << 1];
                return pi[0] == pj[0] && pi[1] == pj[1];
            });
        sorted_pt_inds.erase(it, sorted_pt_inds.end());
        n_points = static_cast<uint32_t>(sorted_pt_inds.size());
    }
    // resort the array of points to accommodate alternating cuts
    alternate_axes(points, &sorted_pt_inds[0], n_points, 0);

    Mesh mesh{points, {}};
    HEdge hull_left, hull_right;
    div_conq_recurse(&mesh, sorted_pt_inds, 0, 0, n_points, hull_left, hull_right);
    std::vector<bool> ghost(mesh.triangles.size(), false);
    *n_triangles = static_cast<uint32_t>(mesh.triangles.size()) - mark_ghost(mesh.triangles, hull_left, ghost);
    form_skeleton(&mesh, n_points, ghost, segments, n_segments);

    auto triangle_indices = std::make_unique<uint32_t[]>(*n_triangles * 3);
    for (uint32_t i = 0, j = 0; i < mesh.triangles.size(); i++) {
        if (ghost[i]) {
            continue;
        }
        const auto& tri = mesh.triangles[i];
        auto data = &triangle_indices[j * 3];
        data[0] = tri.data[0];
        data[1] = tri.data[1];
        data[2] = tri.data[2];
        j += 1;
    }
    return triangle_indices.release();
}
#ifdef EMSCRIPTEN
}
#endif
