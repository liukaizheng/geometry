#include <polygonalization/conforming_mesh.h>
#include <predicates/generic_point.h>

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
    UNDEFINED = 0,
    INTERSECTION = 1,
    IMPROPER_INTERSECTION = 2,
    IMPROPER_INTERSECTION_COUNTED = 3,
    PROPER_INTERSECTION_COUNTED = 4,
    OVERLAP2D_F0 = 10,
    OVERLAP2D_F1 = 11,
    OVERLAP2D_F2 = 12,
    OVERLAP2D_F3 = 13,
    OVERLAP2D_F0_COUNTED = 20,
    OVERLAP2D_F1_COUNTED = 21,
    OVERLAP2D_F2_COUNTED = 22,
    OVERLAP2D_F3_COUNTED = 23
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

inline uint32_t tet_vert_on_constraint_sides(
    TetMesh& mesh, const uint32_t va, const uint32_t vb, uint32_t* connect_verts, std::vector<uint32_t>& temp_tets,
    std::vector<uint32_t>& intersected_tets, std::vector<IType>& intersect_marks
) {
    const uint32_t num_inc_tets = mesh.incident(va, temp_tets);
    const auto fill = [connect_verts](auto... vals) {
        int idx = 0;
        for (const auto val : {vals...}) {
            connect_verts[idx++] = val;
        }
    };

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
            fill(static_cast<uint32_t>(1), vb);
            break;
        }
        // the vertices of opposite face aganist "va"
        std::array<uint32_t, 3> oppo_verts;
        for (uint32_t j = 0, k = 0; j < 4; j++) {
            if (tet.data[j] != va) {
                oppo_verts[k++] = tet.data[j];
            }
        }

        // segemnt and triangle properly intersect
        if (vert_inner_segment_cross_inner_triangle(mesh, va, vb, oppo_verts[0], oppo_verts[1], oppo_verts[2])) {
            fill(static_cast<uint32_t>(3), oppo_verts[0], oppo_verts[1], oppo_verts[2]);
            break;
        }

        // segment and segment properly intersect
        for (uint32_t j = 0; j < 3; j++) {
            const uint32_t ua = oppo_verts[j];
            const uint32_t ub = oppo_verts[(j + 1) % 3];
            if (vert_inner_segments_cross(mesh, va, vb, ua, ub)) {
                fill(static_cast<uint32_t>(2), ua, ub);
                break;
            }
        }

        // tet vertex on segment
        for (uint32_t j = 0; j < 3; j++) {
            if (vert_point_in_inner_segment(mesh, oppo_verts[j], va, vb)) {
                fill(static_cast<uint32_t>(1), oppo_verts[j]);
                break;
            }
        }
    }
    return tid;
}

inline bool verts_in_same_half_space(
    const TetMesh& mesh, const uint32_t u1, const uint32_t u2, const uint32_t v1, const uint32_t v2, const uint32_t v3
) {
    return GenericPoint3D::sign_orient3d(mesh.point(v1), mesh.point(v2), mesh.point(v3), mesh.point(u1)) ==
           GenericPoint3D::sign_orient3d(mesh.point(v1), mesh.point(v2), mesh.point(v3), mesh.point(u2));
}

inline void tet_edge_cross_constraint_side(
    const TetMesh& mesh, const TriFace& edge, const uint32_t va, const uint32_t vb, uint32_t* connect_verts,
    std::vector<uint32_t>& intersected_tets, std::vector<IType>& intersect_marks
) {
    TriFace spin{edge};
    bool found = false;
    do {
        const uint32_t tid = edge.tet;
        if (!mesh.is_hull_tet(tid)) {
            if (intersect_marks[tid] == IType::UNDEFINED) {
                intersected_tets.emplace_back(tid);
                intersect_marks[tid] = IType::INTERSECTION;
            }

            if (!found) {
                const uint32_t oppo_verts[2]{mesh.apex(spin), mesh.oppo(spin)};
                if (oppo_verts[0] == vb || oppo_verts[1] == vb) {
                    mesh.fnext_self(spin);
                    found = true;
                    continue;
                }
            }
        }
        mesh.fnext_self(spin);
    } while (spin.tet != edge.tet);
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
            mesh, triangle[va], triangle[vb], connect_verts.data(), temp_tets, intersected_tets, intersect_mark
        );
        while (connect_verts[1] != vb) {
            switch (connect_verts[0]) {
            case 1:
                tid = tet_vert_on_constraint_sides(
                    mesh, connect_verts[1], vb, connect_verts.data(), temp_tets, intersected_tets, intersect_mark
                );
                break;
            case 2:
                TriFace edge(tid, TetMesh::vpivot[mesh.tets[tid].index(connect_verts[1])]);
                // the third must be succeed
                for (uint32_t j = 0; j < 2; j++) {
                    if (mesh.dest(edge) == connect_verts[2]) {
                        break;
                    }
                    edge.eprev_esym_self();
                }
                tet_edge_cross_constraint_side(
                    mesh, edge, va, vb, connect_verts.data(), intersected_tets, intersect_mark
                );
                break;
            }
        }
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
    }
}
