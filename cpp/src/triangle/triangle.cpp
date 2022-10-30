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

// extern "C" {
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

    const auto iter = std::remove_if(mesh.triangles.begin(), mesh.triangles.end(), [](const Triangle& tri) {
        return tri.data[0] == INVALID || tri.data[1] == INVALID || tri.data[2] == INVALID;
    });
    *n_triangles = static_cast<uint32_t>(std::distance(mesh.triangles.begin(), iter));
    auto triangle_indices = std::make_unique<uint32_t[]>(*n_triangles * 3);
    for (uint32_t i = 0; i < *n_triangles; i++) {
        auto data = &triangle_indices[i * 3];
        data[0] = mesh.triangles[i].data[0];
        data[1] = mesh.triangles[i].data[1];
        data[2] = mesh.triangles[i].data[2];
    }
    return triangle_indices.release();
}
//}
