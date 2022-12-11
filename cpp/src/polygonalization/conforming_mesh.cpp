#include <polygonalization/conforming_mesh.h>
#include <predicates/generic_point.h>

#include <algorithm>
#include <array>
#include <unordered_map>
#include <vector>

using Edges = std::vector<std::unordered_map<uint32_t, std::vector<std::array<uint32_t, 2>>>>;

inline Edges make_edges(const Constraints& constraints, const uint32_t n_points) {
    Edges edges(n_points - 1);
    for (uint32_t i = 0; i < constraints.n_triangles; i++) {
        const uint32_t* triangle = &constraints.triangles[i * 3];
        for (uint32_t j = 0; j < 3; j++) {
            uint32_t pa = triangle[(j + 1) % 3];
            uint32_t pb = triangle[(j + 2) % 3];
            if (pa > pb) {
                std::swap(pa, pb);
            }
            auto& map = edges[pa];
            auto it = map.find(pb);
            if (it != map.end()) {
                it->second.emplace_back(std::array<uint32_t, 2>{{i, j}});
            } else {
                map.emplace(pb, std::vector{std::array<uint32_t, 2>{{i, j}}});
            }
        }
    }
    return edges;
}

inline void add_virtual_constraint(
    const TetMesh& mesh, Constraints& constraints, const std::vector<std::array<uint32_t, 2>>& hedges
) {
    // At least one half edge has been guaranteed
    const auto& first_he = hedges[0];
    const uint32_t* const triangle = &constraints.triangles[first_he[0] * 3];
    uint32_t apex = triangle[first_he[1]];
    for (uint32_t i = 1; i < hedges.size(); i++) {
        const auto& he = hedges[i];
        const uint32_t* tri = &constraints.triangles[he[0] * 3];
        // unplannar: return
        if (mesh.orient3d(apex, tri[0], tri[1], tri[2]) != 0) {
            return;
        }
    }
    if (hedges.size() > 1) {
        const int dom = GenericPoint3D::max_component_at_triangle_normal(
            mesh.point(triangle[0]), mesh.point(triangle[1]), mesh.point(triangle[2])
        );
        auto assign = [dom](const double* src, double* tar) {
            for (int i = 0; i < 2; i++) {
                const int idx = (dom + i + 1) % 3;
                tar[i] = src[idx];
            }
        };
        double pa[2], pb[2], pc[2];
        assign(mesh.point(triangle[(first_he[1] + 1) % 3]), pa);
        assign(mesh.point(triangle[(first_he[1] + 2) % 3]), pb);
        assign(mesh.point(apex), pc);
        const double base_ori = orient2d(pa, pb, pc);
        const int base_sign = (base_ori > 0.0) - (base_ori < 0.0);
        for (uint32_t i = 1; i < hedges.size(); i++) {
            const auto& he = hedges[i];
            apex = constraints.triangles[he[0] * 3 + he[1]];
            assign(mesh.point(apex), pc);
            const double new_ori = orient2d(pa, pb, pc);
            const int new_sign = (new_ori > 0.0) - (new_ori < 0.0);
            if (new_sign != base_sign) {
                return;
            }
        }
    }

    // all triangles are on the same side and coplannar
    const uint32_t tid = mesh.p2t[triangle[0]];
    const uint32_t* const tet = mesh.tets[tid].data.data();
    apex = tet[0];
    // if the first three times failed, then the next time must succeed
    for (int i = 1; i < 4; i++) {
        if (apex == triangle[0] || apex == triangle[1] || apex == triangle[2] ||
            mesh.orient3d(apex, triangle[0], triangle[1], triangle[2]) == 0) {
            apex = tet[i];
        } else {
            break;
        }
    }
    // 'emplace_back' maybe make vector reallocate, the pointer of 'triangle' points to nothing.
    // so we must copy its value firstly.
    const uint32_t org = triangle[(first_he[1] + 1) % 3];
    const uint32_t dest = triangle[(first_he[1] + 2) % 3];
    constraints.triangles.emplace_back(org);
    constraints.triangles.emplace_back(dest);
    constraints.triangles.emplace_back(apex);
}

void place_virtual_constraints(const TetMesh& mesh, Constraints& constraints) {
    auto edges = make_edges(constraints, mesh.n_points);
    for (auto&& map : std::move(edges)) {
        if (map.empty()) {
            continue;
        }
        for (auto&& pair : std::move(map)) {
            add_virtual_constraint(mesh, constraints, pair.second);
        }
    }
}

inline TriFace triangle_at_tet(TetMesh& mesh, const uint32_t* tri, std::vector<uint32_t>& temp_tets) {
    uint32_t count = 0;

    const auto push = [&count, &mesh, &temp_tets](const uint32_t t) {
        if (count < temp_tets.size()) {
            temp_tets[count] = t;
        } else {
            temp_tets.emplace_back(t);
        }
        mesh.mark_test(t);
        count += 1;
    };
    push(mesh.p2t[tri[0]]);

    // find the second point in triangle
    uint32_t vc = 3;
    TriFace edge;
    for (uint32_t i = 0; i < count; i++) {
        const uint32_t tid = temp_tets[i];
        const uint32_t vid = mesh.tets[tid].index(tri[0]);
        TriFace cur(tid, TetMesh::vpivot[vid]);
        bool stop = false;
        for (uint32_t j = 0; j < 3; j++) {
            const uint32_t dest = mesh.dest(cur);
            if (dest == tri[1] || dest == tri[2]) {
                edge = cur;
                vc = dest == tri[1] ? 2 : 1;
                stop = true;
                break;
            }
            cur.eprev_esym_self();
        }
        if (stop) break;
        for (uint32_t j = 0; j < 4; j++) {
            const uint32_t nei = mesh.tets[tid].nei[j].tet;
            if (mesh.is_hull_tet(nei) || mesh.mark_tested(nei)) {
                continue;
            }
            push(nei);
        }
    }

    for (uint32_t i = 0; i < count; i++) {
        mesh.unmark_test(temp_tets[i]);
    }

    if (vc == 3) {
        return {};
    }

    TriFace spin(edge);
    do {
        if (mesh.apex(spin) == tri[vc]) {
            return spin;
        }
        mesh.fnext_self(spin);
    } while (spin.tet != edge.tet);
    return {};
}

// Intersection Type
enum class IType {
    OVERLAP2D_F0 = 0,
    OVERLAP2D_F1 = 1,
    OVERLAP2D_F2 = 2,
    OVERLAP2D_F3 = 3,
    IMPROPER_INTERSECTION = 4,
    INTERSECTION = 5,
    UNDEFINED = 6
};

inline bool vert_inner_segment_cross_inner_triangle(
    const TetMesh& mesh, const uint32_t u1, const uint32_t u2, const uint32_t v1, const uint32_t v2, const uint32_t v3
) {
    if (u1 == v1 || u1 == v2 || u1 == v3 || u2 == v1 || u2 == v2 || u2 == v3) {
        return false;
    }
    return GenericPoint3D::inner_segment_cross_inner_triangle(
        mesh.point(u1), mesh.point(u2), mesh.point(v1), mesh.point(v2), mesh.point(v3)
    );
}

inline bool
vert_inner_segment_cross_triangle(const TetMesh& mesh, const uint32_t u1, const uint32_t u2, const uint32_t* tri) {
    if (u1 == tri[0] || u1 == tri[1] || u1 == tri[2] || u2 == tri[0] || u2 == tri[1] || u2 == tri[2]) {
        return false;
    }
    return GenericPoint3D::inner_segment_cross_triangle(
        mesh.point(u1), mesh.point(u2), mesh.point(tri[0]), mesh.point(tri[1]), mesh.point(tri[2])
    );
}

inline bool vert_inner_segments_cross(
    const TetMesh& mesh, const uint32_t u1, const uint32_t u2, const uint32_t v1, const uint32_t v2
) {
    if (u1 == v1 || u1 == v2 || u2 == v1 || u2 == v2) {
        return false;
    }
    return GenericPoint3D::inner_segments_cross(mesh.point(u1), mesh.point(u2), mesh.point(v1), mesh.point(v2));
}

inline bool vert_point_in_inner_segment(const TetMesh& mesh, const uint32_t u, const uint32_t v1, const uint32_t v2) {
    return u != v1 && u != v2 && GenericPoint3D::point_in_inner_segment(mesh.point(u), mesh.point(v1), mesh.point(v2));
}

inline bool vert_point_in_segment(const TetMesh& mesh, const uint32_t u, const uint32_t v1, const uint32_t v2) {
    return u == v1 || u == v2 || GenericPoint3D::point_in_inner_segment(mesh.point(u), mesh.point(v1), mesh.point(v2));
}

inline bool vert_point_in_inner_triangle(const TetMesh& mesh, const uint32_t p, const uint32_t* tri) {
    return p != tri[0] && p != tri[1] && p != tri[2] &&
           GenericPoint3D::point_in_inner_triangle(
               mesh.point(p), mesh.point(tri[0]), mesh.point(tri[1]), mesh.point(tri[2])
           );
}

inline bool vert_point_in_triangle(const TetMesh& mesh, const uint32_t p, const uint32_t* tri) {
    return p == tri[0] || p == tri[1] || p == tri[2] ||
           GenericPoint3D::point_in_triangle(mesh.point(p), mesh.point(tri[0]), mesh.point(tri[1]), mesh.point(tri[2]));
}

template <typename... T>
void fill_arr(uint32_t* arr, T... vals) {
    int idx = 0;
    for (const auto val : {vals...}) {
        arr[idx++] = val;
    }
}

inline uint32_t tet_vert_on_constraint_sides(
    TetMesh& mesh, const uint32_t va, const uint32_t vb, uint32_t* connect_verts, std::vector<uint32_t>& temp_tets,
    std::vector<uint32_t>& intersected_tets, std::vector<IType>& intersect_marks
) {
    const uint32_t num_inc_tets = mesh.incident(va, temp_tets);

    // mark intersection
    for (uint32_t i = 0; i < num_inc_tets; i++) {
        const uint32_t t = temp_tets[i];
        if (intersect_marks[t] == IType::UNDEFINED) {
            intersected_tets.emplace_back(t);
            intersect_marks[t] = IType::INTERSECTION;
        }
    }

    // find intersections
    uint32_t tid = TriFace::INVALID;
    for (uint32_t i = 0; i < num_inc_tets; i++) {
        tid = temp_tets[i];
        const Tet& tet = mesh.tets[tid];
        // segment is one of edges of tet
        if (tet.index(vb) != 4) {
            fill_arr(connect_verts, static_cast<uint32_t>(1), vb);
            return tid;
        }
        // the vertices of opposite face aganist "va"
        std::array<uint32_t, 3> oppo_verts;
        for (uint32_t j = 0, k = 0; j < 4; j++) {
            if (tet.data[j] != va) {
                oppo_verts[k++] = tet.data[j];
            }
        }

        // segment and triangle properly intersect
        if (vert_inner_segment_cross_inner_triangle(mesh, va, vb, oppo_verts[0], oppo_verts[1], oppo_verts[2])) {
            fill_arr(connect_verts, static_cast<uint32_t>(3), oppo_verts[0], oppo_verts[1], oppo_verts[2]);
            return tid;
        }

        // segment and segment properly intersect
        for (uint32_t j = 0; j < 3; j++) {
            const uint32_t ua = oppo_verts[j];
            const uint32_t ub = oppo_verts[(j + 1) % 3];
            if (vert_inner_segments_cross(mesh, va, vb, ua, ub)) {
                fill_arr(connect_verts, static_cast<uint32_t>(2), ua, ub);
                return tid;
            }
        }

        // tet vertex on segment
        for (uint32_t j = 0; j < 3; j++) {
            if (vert_point_in_inner_segment(mesh, oppo_verts[j], va, vb)) {
                fill_arr(connect_verts, static_cast<uint32_t>(1), oppo_verts[j]);
                return tid;
            }
        }
    }
    return tid;
}

inline bool verts_in_same_half_space(
    const TetMesh& mesh, const uint32_t u1, const uint32_t u2, const uint32_t v1, const uint32_t v2, const uint32_t v3
) {
    return mesh.orient3d(u1, v1, v2, v3) != mesh.orient3d(u2, v1, v2, v3);
}

inline bool
verts_same_half_plane(const TetMesh& mesh, const uint32_t u1, const uint32_t u2, const uint32_t v1, const uint32_t v2) {
    if (u1 == v1 || u1 == v2 || u2 == v1 || u2 == v2) {
        return false;
    }
    return GenericPoint3D::same_half_plane(mesh.point(u1), mesh.point(u2), mesh.point(v1), mesh.point(v2));
}

inline uint32_t tet_edge_cross_constraint_side(
    const TetMesh& mesh, const TriFace& edge, const uint32_t ea, const uint32_t eb, uint32_t* connect_verts,
    std::vector<uint32_t>& intersected_tets, std::vector<IType>& intersect_marks
) {
    TriFace spin{edge};
    mesh.fnext_self(spin);
    bool found = false;
    const uint32_t va = mesh.org(edge);
    const uint32_t vb = mesh.dest(edge);
    uint32_t next_tet{TriFace::INVALID};
    const auto post_process = [&mesh, &spin, &found, &next_tet](const uint32_t t) {
        mesh.fnext_self(spin);
        next_tet = t;
        found = true;
    };
    while (spin.tet != edge.tet) {
        const uint32_t tid = spin.tet;
        if (!mesh.is_hull_tet(tid)) {
            if (intersect_marks[tid] == IType::UNDEFINED) {
                intersected_tets.emplace_back(tid);
                intersect_marks[tid] = IType::INTERSECTION;
            }

            if (!found) {
                const uint32_t vc = mesh.apex(spin);
                const uint32_t vd = mesh.oppo(spin);
                // segment is one of edges of tet
                if (vc == eb || vd == eb) {
                    fill_arr(connect_verts, static_cast<uint32_t>(1), eb);
                    post_process(tid);
                    continue;
                }
                // segment and triangle properly intersect
                if (vert_inner_segment_cross_inner_triangle(mesh, ea, eb, vc, vd, va) &&
                    verts_in_same_half_space(mesh, vb, ea, vc, vd, va)) {
                    fill_arr(connect_verts, static_cast<uint32_t>(3), vc, vd, va);
                    post_process(tid);
                    continue;
                }
                if (vert_inner_segment_cross_inner_triangle(mesh, ea, eb, vc, vd, vb) &&
                    verts_in_same_half_space(mesh, va, ea, vc, vd, vb)) {
                    fill_arr(connect_verts, static_cast<uint32_t>(3), vc, vd, vb);
                    post_process(tid);
                    continue;
                }
                // segment and segment properly intersect
                if (vert_inner_segments_cross(mesh, ea, eb, vc, vd) &&
                    verts_in_same_half_space(mesh, ea, va, vc, vd, vb)) {
                    fill_arr(connect_verts, static_cast<uint32_t>(2), vc, vd);
                    post_process(tid);
                    continue;
                }
                // tet vertex on segment
                if (vert_point_in_inner_segment(mesh, vc, ea, eb) && verts_same_half_plane(mesh, vc, eb, va, vb)) {
                    fill_arr(connect_verts, static_cast<uint32_t>(1), vc);
                    post_process(tid);
                    continue;
                }
                if (vert_point_in_inner_segment(mesh, vd, ea, eb) && verts_same_half_plane(mesh, vd, eb, va, vb)) {
                    fill_arr(connect_verts, static_cast<uint32_t>(1), vd);
                    post_process(tid);
                    continue;
                }

                // segment and segment properly intersect
                for (const uint32_t v1 : {va, vb}) {
                    for (const uint32_t v2 : {vc, vd}) {
                        if (vert_inner_segments_cross(mesh, ea, eb, v1, v2) &&
                            verts_same_half_plane(mesh, v2, eb, va, vb)) {
                            fill_arr(connect_verts, static_cast<uint32_t>(2), v1, v2);
                            post_process(tid);
                            break;
                        }
                    }
                    if (found) {
                        continue;
                    }
                }
            }
        }
        mesh.fnext_self(spin);
    }
    return next_tet;
}

inline uint32_t tet_face_pierced_constraint_side(
    const TetMesh& mesh, TriFace&& tri, const uint32_t ea, uint32_t eb, uint32_t* connect_verts,
    std::vector<uint32_t>& intersected_tets, std::vector<IType>& intersect_marks
) {
    const TriFace& next_tet = mesh.tets[tri.tet].nei[tri.ver & 3];
    const uint32_t next_tid = next_tet.tet;
    if (intersect_marks[next_tid] == IType::UNDEFINED) {
        intersect_marks[next_tid] = IType::INTERSECTION;
        intersected_tets.emplace_back(next_tid);
    }
    const uint32_t v_oppo = mesh.oppo(next_tet);
    if (v_oppo == eb) {
        fill_arr(connect_verts, static_cast<uint32_t>(1), eb);
        return next_tid;
    }

    for (int _ = 0; _ < 3; _++) {
        const uint32_t va = mesh.org(tri);
        const uint32_t vb = mesh.dest(tri);
        if (vert_inner_segment_cross_inner_triangle(mesh, ea, eb, va, vb, v_oppo)) {
            fill_arr(connect_verts, static_cast<uint32_t>(3), va, vb, v_oppo);
            return next_tid;
        }
        tri.enext_self();
    }

    for (const uint32_t v : {mesh.org(tri), mesh.dest(tri), mesh.apex(tri)}) {
        if (vert_inner_segments_cross(mesh, ea, eb, v, v_oppo)) {
            fill_arr(connect_verts, static_cast<uint32_t>(2), v, v_oppo);
            return next_tid;
        }
    }
    // it must be that [ea, eb] pass through v_oppo
    fill_arr(connect_verts, static_cast<uint32_t>(1), v_oppo);
    return next_tid;
}

inline void intersection_constraint_sides(
    TetMesh& mesh, const uint32_t* triangle, std::vector<uint32_t>& temp_tets, std::vector<uint32_t>& intersected_tets,
    std::vector<IType>& intersect_mark
) {
    std::array<uint32_t, 4> connect_verts; // [num, ...verts]
    for (uint32_t i = 0; i < 3; i++) {
        const uint32_t va = triangle[(i + 1) % 3];
        const uint32_t vb = triangle[(i + 2) % 3];
        uint32_t tid = tet_vert_on_constraint_sides(
            mesh, va, vb, connect_verts.data(), temp_tets, intersected_tets, intersect_mark
        );
        while (connect_verts[1] != vb) {
            switch (connect_verts[0]) {
            case 1:
                tid = tet_vert_on_constraint_sides(
                    mesh, connect_verts[1], vb, connect_verts.data(), temp_tets, intersected_tets, intersect_mark
                );
                break;
            case 2: {
                TriFace edge(tid, TetMesh::vpivot[mesh.tets[tid].index(connect_verts[1])]);
                // the third must be succeed
                for (uint32_t j = 0; j < 2; j++) {
                    if (mesh.dest(edge) == connect_verts[2]) {
                        break;
                    }
                    edge.eprev_esym_self();
                }
                tid = tet_edge_cross_constraint_side(
                    mesh, edge, va, vb, connect_verts.data(), intersected_tets, intersect_mark
                );
                break;
            }
            case 3: {
                TriFace edge(tid, TetMesh::vpivot[mesh.tets[tid].index(connect_verts[1])]);
                for (uint32_t j = 0; j < 2; j++) {
                    const uint32_t dest = mesh.dest(edge);
                    const uint32_t apex = mesh.apex(edge);
                    if ((dest == connect_verts[2] && apex == connect_verts[3]) ||
                        (dest == connect_verts[3] && apex == connect_verts[2])) {
                        break;
                    }
                    edge.eprev_esym_self();
                }
                tid = tet_face_pierced_constraint_side(
                    mesh, std::move(edge), va, vb, connect_verts.data(), intersected_tets, intersect_mark
                );
            } break;
            }
        }
    }
}

inline bool intersection_class_tetface(const TetMesh& mesh, const uint32_t* c, const uint32_t* t) {
    if (vert_point_in_triangle(mesh, t[0], c) && vert_point_in_triangle(mesh, t[1], c) &&
        vert_point_in_triangle(mesh, t[2], c)) {
        return true;
    }
    if (vert_inner_segments_cross(mesh, t[0], t[1], c[0], c[1])) return true;
    if (vert_inner_segments_cross(mesh, t[0], t[1], c[1], c[2])) return true;
    if (vert_inner_segments_cross(mesh, t[0], t[1], c[2], c[0])) return true;

    if (vert_inner_segments_cross(mesh, t[1], t[2], c[0], c[1])) return true;
    if (vert_inner_segments_cross(mesh, t[1], t[2], c[1], c[2])) return true;
    if (vert_inner_segments_cross(mesh, t[1], t[2], c[2], c[0])) return true;

    if (vert_inner_segments_cross(mesh, t[2], t[0], c[0], c[1])) return true;
    if (vert_inner_segments_cross(mesh, t[2], t[0], c[1], c[2])) return true;
    if (vert_inner_segments_cross(mesh, t[2], t[0], c[2], c[0])) return true;
    return false;
}

inline bool intersection_class_tetedge(
    const TetMesh& mesh, const uint32_t* c, const uint32_t* coplanar_v, const uint32_t* unplanar_v
) {
    bool va_in = vert_point_in_inner_triangle(mesh, coplanar_v[0], c);
    bool vb_in = vert_point_in_inner_triangle(mesh, coplanar_v[1], c);
    // va and vb both in the interior of triangle, intersection must be improper.
    if (va_in && vb_in) {
        return true;
    }
    int va_on_edge = static_cast<int>(vert_point_in_segment(mesh, coplanar_v[0], c[1], c[2])) |
                     (vert_point_in_segment(mesh, coplanar_v[0], c[2], c[0]) << 1) |
                     (vert_point_in_segment(mesh, coplanar_v[0], c[0], c[1]) << 2);
    int vb_on_edge = static_cast<int>(vert_point_in_segment(mesh, coplanar_v[1], c[1], c[2])) |
                     (vert_point_in_segment(mesh, coplanar_v[1], c[2], c[0]) << 1) |
                     (vert_point_in_segment(mesh, coplanar_v[2], c[0], c[1]) << 2);
    va_in |= va_on_edge > 0;
    vb_in |= vb_on_edge > 0;
    if (va_in && vb_in) {
        // va and vb are not on same edge
        const int on_same_edge = va_on_edge & vb_on_edge;
        if (on_same_edge == 0) {
            return true;
        }
        // edge [vc vd] pierces triangle
        if (vert_inner_segment_cross_triangle(mesh, unplanar_v[0], unplanar_v[1], c)) {
            return true;
        }
        const int bit = on_same_edge == 1 ? 0 : (on_same_edge == 2 ? 1 : 2);
        const uint32_t ea = c[(bit + 1) % 3];
        const uint32_t eb = c[(bit + 2) % 3];
        for (const uint32_t ev : {ea, eb}) {
            for (const uint32_t tv : {coplanar_v[0], coplanar_v[1]}) {
                if (vert_inner_segment_cross_inner_triangle(mesh, ev, c[bit], unplanar_v[0], unplanar_v[1], tv)) {
                    return true;
                }
            }
        }
        // otherwise.. the intersection is proper
        return false;
    } else if (va_in && vb_in) {
        // edge [va vb] crosses triangle edge
        if (vert_inner_segments_cross(mesh, coplanar_v[0], coplanar_v[1], c[0], c[1]) ||
            vert_inner_segments_cross(mesh, coplanar_v[0], coplanar_v[1], c[1], c[2]) ||
            vert_inner_segments_cross(mesh, coplanar_v[0], coplanar_v[1], c[2], c[0])) {
            return true;
        }
        // edge [vc vd] pierces triangle
        if (vert_inner_segment_cross_triangle(mesh, unplanar_v[0], unplanar_v[1], c)) {
            return true;
        }
        for (int i = 0; i < 3; i++) {
            const uint32_t ea = c[i];
            const uint32_t eb = c[(i + 1) % 3];
            for (const uint32_t tv : {coplanar_v[0], coplanar_v[1]}) {
                if (vert_inner_segment_cross_inner_triangle(mesh, ea, eb, unplanar_v[0], unplanar_v[1], tv)) {
                    return true;
                }
            }
        }
        // otherwise.. the intersection is proper
        return false;
    } else {
        // otherwise.. since an intersection exists, it must be improper
        return true;
    }
}

inline bool intersection_class_tetvert(
    const TetMesh& mesh, const uint32_t* c, const uint32_t tv, const uint32_t* under_v, const uint32_t over_v
) {
    if (vert_point_in_inner_triangle(mesh, tv, c)) {
        return true;
    }
    // clang-format off
    bool vert_on_segments = vert_point_in_segment(mesh, tv, c[0], c[1]) ||
                            vert_point_in_segment(mesh, tv, c[1], c[2]) ||
                            vert_point_in_segment(mesh, tv, c[2], c[0]);
    // clang-format on
    // tv isn't on triangle, but intersection exist, sot it's improper
    if (!vert_on_segments) {
        return true;
    }
    // now tv is on a boundary of triangle
    if (vert_inner_segment_cross_triangle(mesh, over_v, under_v[0], c) ||
        vert_inner_segment_cross_triangle(mesh, over_v, under_v[1], c)) {
        return true;
    }
    for (uint32_t i = 0; i < 3; i++) {
        const uint32_t ca = c[i];
        const uint32_t cb = c[(i + 1) % 3];
        for (const uint32_t v : {under_v[0], under_v[1]}) {
            if (vert_inner_segment_cross_inner_triangle(mesh, ca, cb, tv, over_v, v)) {
                return true;
            }
        }
        if (vert_inner_segment_cross_inner_triangle(mesh, ca, cb, tv, over_v, under_v[0]) ||
            vert_inner_segment_cross_inner_triangle(mesh, ca, cb, tv, over_v, under_v[1]) ||
            vert_inner_segment_cross_inner_triangle(mesh, ca, cb, over_v, under_v[0], under_v[1])) {
            return true;
        }
    }
    // otherwise... proper intersection
    return false;
}


inline void find_improper_intersection(
    const TetMesh& mesh, const uint32_t* triangle, const std::vector<uint32_t>& tets, IType* tet_marks
) {
    std::unordered_map<uint32_t, int> vert_ori_map;
    const auto vert_ori = [&vert_ori_map, &triangle, &mesh](const uint32_t vid) {
        const auto iter = vert_ori_map.find(vid);
        if (iter != vert_ori_map.end()) {
            return iter->second;
        } else {
            const int ori = mesh.orient3d(triangle[0], triangle[1], triangle[2], vid);
            vert_ori_map.emplace(vid, ori);
            return ori;
        }
    };

    std::vector<uint32_t> zero; // coplanar vertices
    std::vector<uint32_t> neg;  // vertices with negative orientation
    std::vector<uint32_t> pos;  // vertices with positive orientation
    for (const uint32_t tid : tets) {
        const Tet& tet = mesh.tets[tid];
        int tet_vert_oris[4]{
            vert_ori(tet.data[0]), vert_ori(tet.data[1]), vert_ori(tet.data[2]), vert_ori(tet.data[3])};
        zero.clear();
        neg.clear();
        pos.clear();
        for (uint32_t i = 0; i < 4; i++) {
            if (tet_vert_oris[i] == 0) {
                zero.emplace_back(tet.data[i]);
            } else if (tet_vert_oris[i] < 0) {
                neg.emplace_back(tet.data[i]);
            } else {
                pos.emplace_back(tet.data[i]);
            }
        }
        if (zero.size() == 3) {
            const uint32_t oppo_vid = tet.index(neg.empty() ? pos[0] : neg[0]);
            if (intersection_class_tetface(mesh, triangle, zero.data())) {
                tet_marks[tid] = static_cast<IType>(oppo_vid);
                continue;
            }
        } else if (zero.size() == 2) {
            if (neg.size() != 1) {
                continue;
            }
            const uint32_t unplanar_v[2]{neg[0], pos[0]};
            if (intersection_class_tetedge(mesh, triangle, zero.data(), unplanar_v)) {
                tet_marks[tid] = IType::IMPROPER_INTERSECTION;
                continue;
            }
        } else if (zero.size() == 1) {
            if (neg.empty() || pos.empty()) {
                continue;
            }
            const uint32_t* under_v;
            uint32_t over_v;
            if (neg.size() > pos.size()) {
                under_v = neg.data();
                over_v = pos[0];
            } else {
                under_v = pos.data();
                over_v = neg[0];
            }
            if (intersection_class_tetvert(mesh, triangle, zero[0], under_v, over_v)) {
                tet_marks[tid] = IType::IMPROPER_INTERSECTION;
                continue;
            }
        }
        // since intersection exists, it's improper
        tet_marks[tid] = IType::IMPROPER_INTERSECTION;
    }
}

inline IType tet_intersects_triangle_interior(const TetMesh& mesh, const uint32_t* c, const uint32_t tid) {
    const uint32_t* verts = mesh.tets[tid].data.data();
    int t_ori[4];
    bool t_in_c[4];
    for (uint32_t i = 0; i < 3; i++) {
        t_ori[i] = mesh.orient3d(verts[i], c[0], c[1], c[2]);
        t_in_c[i] = t_ori[i] == 0 && vert_point_in_inner_triangle(mesh, verts[i], c);
    }
    if (mesh.is_hull_tet(tid)) {
        if (t_in_c[0] && t_in_c[1] && t_in_c[2]) {
            return IType::INTERSECTION;
        } else {
            return IType::UNDEFINED;
        }
    }
    t_ori[3] = mesh.orient3d(verts[3], c[0], c[1], c[2]);
    t_in_c[3] = t_ori[3] == 0 && vert_point_in_inner_triangle(mesh, verts[3], c);
    std::vector<uint32_t> verts_in_c;
    std::vector<uint32_t> verts_out_c;
    for (uint32_t i = 0; i < 4; i++) {
        if (t_in_c[0]) {
            verts_in_c.emplace_back(i);
        } else {
            verts_out_c.emplace_back(i);
        }
    }
    if (verts_in_c.size() == 3) {
        return static_cast<IType>(verts_out_c[0]);
    } else if (verts_in_c.size() == 2) {
        if (t_ori[verts_out_c[0]] == t_ori[verts_out_c[1]]) {
            return IType::INTERSECTION;
        } else {
            return IType::IMPROPER_INTERSECTION;
        }
    } else if (verts_in_c.size() == 1) {
        if (t_ori[verts_out_c[0]] == t_ori[verts_out_c[1]] && t_ori[verts_out_c[0]] == t_ori[verts_out_c[2]]) {
            return IType::INTERSECTION;
        } else {
            return IType::IMPROPER_INTERSECTION;
        }
    } else {
        for (uint32_t i = 0; i < 3; i++) {
            for (uint32_t j = i + 1; j < 4; j++) {
                if (t_ori[i] != t_ori[j] &&
                    vert_inner_segment_cross_inner_triangle(mesh, verts[i], verts[j], c[0], c[1], c[2])) {
                    return IType::IMPROPER_INTERSECTION;
                }
            }
        }
        return IType::UNDEFINED;
    }
}

inline void in_constraint_interior_test(
    TetMesh& mesh, const uint32_t* triangle, const uint32_t tid, std::vector<uint32_t>& intersected_tets,
    std::vector<uint32_t>& no_intersect_tets, IType* marks
) {
    const auto type = tet_intersects_triangle_interior(mesh, triangle, tid);
    if (type == IType::UNDEFINED) {
        no_intersect_tets.emplace_back(tid);
        mesh.mark_test(tid);
    } else {
        marks[tid] = type;
        intersected_tets.emplace_back(tid);
    }
}

inline void intersection_constraint_interior(
    TetMesh& mesh, const uint32_t* triangle, std::vector<uint32_t>& intersected_tets, IType* intersect_marks
) {
    std::vector<uint32_t> no_intersection_tets;
    for (uint32_t i = 0; i < intersected_tets.size(); i++) {
        const uint32_t tid = intersected_tets[i];
        for (uint32_t j = 0; j < 4; j++) {
            const uint32_t adj_tet = mesh.tets[tid].nei[j].tet;
            if (intersect_marks[adj_tet] == IType::UNDEFINED && !mesh.mark_tested(adj_tet)) {
                in_constraint_interior_test(
                    mesh, triangle, adj_tet, intersected_tets, no_intersection_tets, intersect_marks
                );
            }
        }
    }
    for (const uint32_t t : no_intersection_tets) {
        mesh.unmark_test(t);
    }
}

void insert_constraints(TetMesh& mesh, const Constraints& constraints, std::vector<std::vector<uint32_t>>* tet_map) {
    const uint32_t n_triangles = static_cast<uint32_t>(constraints.triangles.size() / 3);
    std::vector<uint32_t> temp_tets;
    std::vector<IType> intersect_marks(mesh.tets.size(), IType::UNDEFINED);
    for (uint32_t i = 0; i < n_triangles; i++) {
        const uint32_t* triangle = &constraints.triangles[i * 3];
        const TriFace tet_face = triangle_at_tet(mesh, triangle, temp_tets);
        if (tet_face.tet != TriFace::INVALID) {
            tet_map[tet_face.ver & 3][tet_face.tet].emplace_back(i);
            const TriFace& nei = mesh.tets[tet_face.tet].nei[tet_face.ver & 3];
            tet_map[nei.ver & 3][nei.tet].emplace_back(i);
            continue;
        }
        std::vector<uint32_t> intersected_tets;
        intersection_constraint_sides(mesh, triangle, temp_tets, intersected_tets, intersect_marks);
        find_improper_intersection(mesh, triangle, intersected_tets, intersect_marks.data());
        intersection_constraint_interior(mesh, triangle, intersected_tets, intersect_marks.data());
        for (const uint32_t tid : intersected_tets) {
            uint32_t type = static_cast<uint32_t>(intersect_marks[tid]);
            if (type < 5) {
                tet_map[type][tid].emplace_back(i);
            }
            intersect_marks[tid] = IType::UNDEFINED;
        }
    }
}
