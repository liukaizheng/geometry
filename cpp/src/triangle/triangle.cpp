#include <memory>
#include <triangle/triangle.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <cstdint>

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
    const std::vector<uint32_t>& sorted_pt_inds;
    std::vector<Triangle> triangles;
};

static void alternate_axes(const double* points, uint32_t* indices, const uint32_t len, uint32_t axis) {
    uint32_t divider = len >> 1;
    if (len <= 3) {
        axis = 0;
    }
    std::nth_element(indices, indices + divider, indices + len, [&points, axis](const uint32_t i, const uint32_t j) {
        const double* pi = &points[i * 2];
        const double* pj = &points[j * 2];
        return pi[axis] < pj[axis];
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
merge_hulls(Mesh* m, const uint32_t axis, HEdge& farleft, HEdge& innerleft, HEdge& innerright, HEdge& farright) {
    uint32_t innerleftdest = dest(m, innerleft);
    uint32_t innerleftapex = apex(m, innerleft);
    uint32_t innerrightorg = org(m, innerright);
    uint32_t innerrightapex = apex(m, innerright);
    uint32_t farleftpt, farleftapex, farrightpt, farrightapex;
    // Special treatment for horizontal cuts
    if (axis == 1) {
        farleftpt = org(m, farleft);
        farleftapex = apex(m, farleft);
        farrightpt = dest(m, farright);
        farrightapex = apex(m, farright);
        // The pointers to the extremal vertices are shifted to point to the
        // topmost and bottommost vertex of each hull, rather than the
        // leftmost and rightmost vertices.
        while (point(m->points, farleftapex)[1] < point(m->points, farleftpt)[1]) {
            next_self(farleft);
            sym_self(m->triangles, farleft);
            farleftpt = farleftapex;
            farleftapex = apex(m, farleft);
        }

        HEdge check_edge;
        sym(m->triangles, innerleft, check_edge);
        uint32_t check_vertex = apex(m, check_edge);
        while (point(m->points, check_vertex)[1] > point(m->points, innerleftdest)[1]) {
            next(check_edge, innerleft);
            innerleftapex = innerleftdest;
            innerleftdest = check_vertex;
            sym(m->triangles, innerleft, check_edge);
            check_vertex = apex(m, check_edge);
        }

        while (point(m->points, innerrightapex)[1] < point(m->points, innerrightorg)[1]) {
            next_self(innerright);
            sym_self(m->triangles, innerright);
            innerrightorg = innerrightapex;
            innerrightapex = apex(m, innerright);
        }

        sym(m->triangles, farright, check_edge);
        check_vertex = apex(m, check_edge);
        while (point(m->points, check_vertex)[1] > point(m->points, farrightpt)[1]) {
            next(check_edge, farright);
            farrightapex = farrightpt;
            farrightpt = check_vertex;
            sym(m->triangles, farright, check_edge);
            check_vertex = apex(m, check_edge);
        }
    }

    // Find a line tangent to and below both hulls
    bool change_made = true;
    while (change_made) {
        change_made = false;
        // Make innerleftdest the "bottommost" vertex of the left hull
        if (counterclockwise(m->points, innerleftdest, innerleftapex, innerrightorg) > 0.0) {
            prev_self(innerleft);
            sym_self(m->triangles, innerleft);
            innerleftdest = innerleftapex;
            innerleftapex = apex(m, innerleft);
            change_made = true;
        }
        // Make innerrightorg the "bottommost" vertex of the right hull
        if (counterclockwise(m->points, innerrightapex, innerrightorg, innerleftdest) > 0.0) {
            next_self(innerright);
            sym_self(m->triangles, innerright);
            innerrightorg = innerrightapex;
            innerrightapex = apex(m, innerright);
            change_made = true;
        }
    }

    // Find the two candidates to be the next "gear tooth"
    HEdge leftcand, rightcand, baseedge;
    sym(m->triangles, innerleft, leftcand);
    sym(m->triangles, innerright, rightcand);
    // Create the bottom new bounding triangle
    make_triangle(m->triangles, baseedge);
    /* Connect it to the bounding boxes of the left and right triangulations. */
    bond(m->triangles, baseedge, innerleft);
    next_self(baseedge);
    bond(m->triangles, baseedge, innerright);
    next_self(baseedge);
    set_org(m->triangles, baseedge, innerrightorg);
    set_dest(m->triangles, baseedge, innerleftdest);


    // Fix the extreme triangles if necessary
    farleftpt = org(m, farleft);
    if (innerleftdest == farleftpt) {
        next(baseedge, farleft);
    }
    farrightpt = dest(m, farright);
    if (innerrightorg == farrightpt) {
        prev(baseedge, farright);
    }

    // The vertices of the current knitting edge
    uint32_t lowerleft = innerleftdest;
    uint32_t lowerright = innerrightorg;
    // The candidate vertices for knitting
    uint32_t upperleft = apex(m, leftcand);
    uint32_t upperright = apex(m, rightcand);
    // Walk up the gap between the two triangulations, knitting them together
    while (true) {
        const auto leftfinished = counterclockwise(m->points, upperleft, lowerleft, lowerright) <= 0.0;
        const auto rightfinished = counterclockwise(m->points, upperright, lowerleft, lowerright) <= 0.0;
        HEdge checkedge, nextedge;
        if (leftfinished && rightfinished) {
            // Create the top new bounding triangle
            make_triangle(m->triangles, nextedge);
            set_org(m->triangles, nextedge, lowerleft);
            set_dest(m->triangles, nextedge, lowerright);
            // Apex is intentionally left INVALID
            // Connect it to the bounding boxes of the two triangulations
            bond(m->triangles, nextedge, baseedge);
            next_self(nextedge);
            bond(m->triangles, nextedge, rightcand);
            next_self(nextedge);
            bond(m->triangles, nextedge, leftcand);

            // Special treatment for horizontal cuts
            if (axis == 1) {
                farleftpt = org(m, farleft);
                farleftapex = apex(m, farleft);
                farrightpt = dest(m, farright);
                farrightapex = apex(m, farright);
                sym(m->triangles, farleft, checkedge);
                uint32_t checkvertex = apex(m, checkedge);
                // The pointers to the extremal vertices are restored to the
                // leftmost and rightmost vertices (rather than topmost and
                // bottommost)
                while (point(m->points, checkvertex)[0] < point(m->points, farleftpt)[0]) {
                    prev(checkedge, farleft);
                    farleftapex = farleftpt;
                    farleftpt = checkvertex;
                    sym(m->triangles, farleft, checkedge);
                    checkvertex = apex(m, checkedge);
                }
                while (point(m->points, farrightapex)[0] > point(m->points, farrightpt)[0]) {
                    prev_self(farright);
                    sym_self(m->triangles, farright);
                    farrightpt = farrightapex;
                    farrightapex = apex(m, farright);
                }
            }
            return;
        }

        // Consider eliminating edges from the left triangulation
        if (!leftfinished) {
            // What vertex would be exposed if an edge were deleted
            prev(leftcand, nextedge);
            sym_self(m->triangles, nextedge);
            uint32_t nextapex = apex(m, nextedge);
            // If nextapex is INVALID, then no vertex would be exposed; the
            // triangulation would have been eaten right through.
            if (nextapex != INVALID) {
                // Check whether the edge is Delaunay
                auto badedge = incircle(m->points, lowerleft, lowerright, upperleft, nextapex) > 0.0;
                HEdge topcasing, sidecasing, outercasing;
                while (badedge) {
                    // Eliminate the edge with an edge flip.  As a result, the
                    // left triangulation will have one more boundary triangle.
                    next_self(nextedge);
                    sym(m->triangles, nextedge, topcasing);
                    next_self(nextedge);
                    sym(m->triangles, nextedge, sidecasing);
                    bond(m->triangles, nextedge, topcasing);
                    bond(m->triangles, leftcand, sidecasing);
                    next_self(leftcand);
                    sym(m->triangles, leftcand, outercasing);
                    prev_self(nextedge);
                    bond(m->triangles, nextedge, outercasing);
                    // Correct the vertices to reflect the edge flip
                    set_org(m->triangles, leftcand, lowerleft);
                    set_dest(m->triangles, leftcand, INVALID);
                    set_apex(m->triangles, leftcand, nextapex);
                    set_org(m->triangles, nextedge, INVALID);
                    set_dest(m->triangles, nextedge, upperleft);
                    set_apex(m->triangles, nextedge, nextapex);
                    // Consider the newly exposed vertex
                    upperleft = nextapex;
                    // What vertex would be exposed if another edge were deleted?
                    copy(sidecasing, nextedge);
                    nextapex = apex(m, nextedge);
                    if (nextapex != INVALID) {
                        // Check whether the edge is Delaunay
                        badedge = incircle(m->points, lowerleft, lowerright, upperleft, nextapex) > 0.0;
                    } else {
                        // Avoid eating right through the triangulation
                        badedge = false;
                    }
                }
            }
        }

        // Consider eliminating edges from the right triangulation
        if (!rightfinished) {
            next(rightcand, nextedge);
            sym_self(m->triangles, nextedge);
            uint32_t nextapex = apex(m, nextedge);
            if (nextapex != INVALID) {
                // Check whether the edge is Delaunay
                bool badedge = incircle(m->points, lowerleft, lowerright, upperright, nextapex) > 0.0;
                HEdge topcasing, sidecasing, outercasing;
                while (badedge) {
                    prev_self(nextedge);
                    sym(m->triangles, nextedge, topcasing);
                    prev_self(nextedge);
                    sym(m->triangles, nextedge, sidecasing);
                    bond(m->triangles, nextedge, topcasing);
                    bond(m->triangles, rightcand, sidecasing);
                    prev_self(rightcand);
                    sym(m->triangles, rightcand, outercasing);
                    next_self(nextedge);
                    bond(m->triangles, nextedge, outercasing);

                    set_org(m->triangles, rightcand, INVALID);
                    set_dest(m->triangles, rightcand, lowerright);
                    set_apex(m->triangles, rightcand, nextapex);
                    set_org(m->triangles, nextedge, upperright);
                    set_dest(m->triangles, nextedge, INVALID);
                    set_apex(m->triangles, nextedge, nextapex);

                    upperright = nextapex;

                    copy(sidecasing, nextedge);
                    nextapex = apex(m, nextedge);
                    if (nextapex != INVALID) {
                        badedge = incircle(m->points, lowerleft, lowerright, upperright, nextapex) > 0.0;
                    } else {
                        badedge = false;
                    }
                }
            }
        }


        if (leftfinished ||
            (!rightfinished && (incircle(m->points, upperleft, lowerleft, lowerright, upperright) > 0.0))) {
            // Knit the triangulations, adding an edge from `lowerleft'
            // to `upperright'
            bond(m->triangles, baseedge, rightcand);
            prev(rightcand, baseedge);
            set_dest(m->triangles, baseedge, lowerleft);
            lowerright = upperright;
            sym(m->triangles, baseedge, rightcand);
            upperright = apex(m, rightcand);
        } else {
            // Knit the triangulations, adding an edge from `upperleft'
            // to `lowerright'
            bond(m->triangles, baseedge, leftcand);
            next(leftcand, baseedge);
            set_org(m->triangles, baseedge, lowerright);
            lowerleft = upperleft;
            sym(m->triangles, baseedge, leftcand);
            upperleft = apex(m, leftcand);
        }
    }
}

static void div_conq_recurse(
    Mesh* m, const uint32_t axis, const uint32_t start, const uint32_t end, HEdge& farleft, HEdge& farright
) {
    const auto len = end - start;
    if (len == 2) {
        make_triangle(m->triangles, farleft);
        set_org(m->triangles, farleft, m->sorted_pt_inds[start]);
        set_dest(m->triangles, farleft, m->sorted_pt_inds[start + 1]);

        make_triangle(m->triangles, farright);
        set_org(m->triangles, farright, m->sorted_pt_inds[start + 1]);
        set_dest(m->triangles, farright, m->sorted_pt_inds[start]);
        bond(m->triangles, farleft, farright);

        prev_self(farleft);
        next_self(farright);
        bond(m->triangles, farleft, farright);

        prev_self(farleft);
        next_self(farright);
        bond(m->triangles, farleft, farright);

        // ensura that the origin of `farleft` is `start`
        prev(farright, farleft);
    } else if (len == 3) {
        HEdge midtri, tri1, tri2, tri3;
        make_triangle(m->triangles, midtri);
        make_triangle(m->triangles, tri1);
        make_triangle(m->triangles, tri2);
        make_triangle(m->triangles, tri3);
        const auto area = counterclockwise(
            m->points, m->sorted_pt_inds[start], m->sorted_pt_inds[start + 1], m->sorted_pt_inds[start + 2]
        );
        if (area == 0.0) {
            // Three collinear vertices; the triangulation is two edges
            // const auto p1 = &m->points[m->sorted_pt_inds[start]];
            // const auto p2 = &m->points[m->sorted_pt_inds[start + 1]];
            // const auto p3 = &m->points[m->sorted_pt_inds[start + 2]];
            set_org(m->triangles, midtri, m->sorted_pt_inds[start]);
            set_dest(m->triangles, midtri, m->sorted_pt_inds[start + 1]);
            set_org(m->triangles, tri1, m->sorted_pt_inds[start + 1]);
            set_dest(m->triangles, tri1, m->sorted_pt_inds[start]);
            set_org(m->triangles, tri2, m->sorted_pt_inds[start + 2]);
            set_dest(m->triangles, tri2, m->sorted_pt_inds[start + 1]);
            set_org(m->triangles, tri3, m->sorted_pt_inds[start + 1]);
            set_dest(m->triangles, tri3, m->sorted_pt_inds[start + 2]);
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

            copy(tri1, farleft);
            copy(tri2, farright);
        } else {

            // The three vertices are not collinear; the triangulation is one
            // triangle, namely `midtri'.
            set_org(m->triangles, midtri, m->sorted_pt_inds[start]);
            set_dest(m->triangles, tri1, m->sorted_pt_inds[start]);
            set_org(m->triangles, tri3, m->sorted_pt_inds[start]);
            // Apices of tri1, tri2, and tri3 are left NULL
            if (area > 0.0) {
                // The vertices are in counterclockwise order
                set_dest(m->triangles, midtri, m->sorted_pt_inds[start + 1]);
                set_org(m->triangles, tri1, m->sorted_pt_inds[start + 1]);
                set_dest(m->triangles, tri2, m->sorted_pt_inds[start + 1]);
                set_apex(m->triangles, midtri, m->sorted_pt_inds[start + 2]);
                set_org(m->triangles, tri2, m->sorted_pt_inds[start + 2]);
                set_dest(m->triangles, tri3, m->sorted_pt_inds[start + 2]);
            } else {
                // The vertices are in clockwise order
                set_dest(m->triangles, midtri, m->sorted_pt_inds[start + 2]);
                set_org(m->triangles, tri1, m->sorted_pt_inds[start + 2]);
                set_dest(m->triangles, tri2, m->sorted_pt_inds[start + 2]);
                set_apex(m->triangles, midtri, m->sorted_pt_inds[start + 1]);
                set_org(m->triangles, tri2, m->sorted_pt_inds[start + 1]);
                set_dest(m->triangles, tri3, m->sorted_pt_inds[start + 1]);
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
            copy(tri1, farleft);
            // Ensure that the destination of `farright' is `start + 2`
            if (area > 0.0) {
                copy(tri2, farright);
            } else {
                next(farleft, farright);
            }
        }
    } else {
        const auto divider = len >> 1;
        HEdge innerleft, innerright;
        div_conq_recurse(m, 1 - axis, start, start + divider, farleft, innerleft);
        div_conq_recurse(m, 1 - axis, start + divider, end, innerright, farright);
        merge_hulls(m, axis, farleft, innerleft, innerright, farright);
    }
}

uint32_t*
triangulate(const double* points, uint32_t n_points, const uint32_t* segments, const uint32_t n_segments, uint32_t* n_triangles) {
    std::vector<uint32_t> sorted_pt_inds(n_points);
    std::iota(sorted_pt_inds.begin(), sorted_pt_inds.end(), 0);
    std::sort(sorted_pt_inds.begin(), sorted_pt_inds.end(), [points](const uint32_t i, const uint32_t j) {
        const double* pi = &points[i * 2];
        const double* pj = &points[j * 2];
        return pi[0] < pj[0] || (pi[0] == pj[0] && pi[1] < pj[1]);
    });
    {
        const auto it =
            std::unique(sorted_pt_inds.begin(), sorted_pt_inds.end(), [points](const uint32_t i, const uint32_t j) {
                const double* pi = &points[i * 2];
                const double* pj = &points[j * 2];
                return pi[0] == pj[0] && pi[1] == pj[1];
            });
        sorted_pt_inds.erase(it, sorted_pt_inds.end());
        n_points = static_cast<uint32_t>(sorted_pt_inds.size());
    }
    // resort the array of points to accommodate alternating cuts
    alternate_axes(points, &sorted_pt_inds[0], n_points, 0);

    Mesh mesh{points, sorted_pt_inds, {}};
    HEdge hullleft, hullright;
    div_conq_recurse(&mesh, 0, 0, n_points, hullleft, hullright);

    const auto iter = std::remove_if(mesh.triangles.begin(), mesh.triangles.end(), [](const Triangle& tri) {
        return tri.data[0] == INVALID || tri.data[1] == INVALID || tri.data[2] == INVALID;
    });
    *n_triangles = static_cast<uint32_t>(std::distance(mesh.triangles.begin(), iter));
    auto triangle_indices = std::make_unique<uint32_t[]>(*n_triangles * 3);
    for(uint32_t i = 0; i < *n_triangles; i++) {
        auto data = &triangle_indices[i * 3];
        data[0] = mesh.triangles[i].data[0];
        data[1] = mesh.triangles[i].data[1];
        data[2] = mesh.triangles[i].data[2];
    }
    return triangle_indices.release();
}

