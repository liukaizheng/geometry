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

static constexpr std::array<uint32_t, 4> vpivot{{11, 8, 9, 10}};

inline TriFace triangle_at_tet(TetMesh& mesh, const uint32_t* tri) {
    uint32_t count = 0;

    const auto push = [&count, &mesh](const uint32_t t) {
        if (count < TetMesh::temp_tets.size()) {
            TetMesh::temp_tets[count] = t;
        } else {
            TetMesh::temp_tets.emplace_back(t);
        }
        mesh.mark_test(t);
        count += 1;
    };
    push(mesh.p2t[tri[0]]);

    // find the second point in triangle
    uint32_t vc = 3;
    TriFace edge;
    for (uint32_t i = 0; i < count; i++) {
        const uint32_t tid = TetMesh::temp_tets[i];
        const uint32_t vid = mesh.tets[tid].index(tri[0]);
        TriFace cur(tid, vpivot[vid]);
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
        mesh.unmark_test(TetMesh::temp_tets[i]);
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

void insert_constraints(
    TetMesh& mesh, const Constraints& constraints, std::vector<std::vector<uint32_t>>* tet_map
) {
    const uint32_t n_triangles = static_cast<uint32_t>(constraints.triangles.size() / 3);
    std::vector<int> tet_marks(mesh.tets.size());
    for (uint32_t i = 0; i < n_triangles; i++) {
        const uint32_t* triangle = &constraints.triangles[i * 3];
        const TriFace tet_face = triangle_at_tet(mesh, triangle);
        if (tet_face.tet != TriFace::INVALID) {
            tet_map[tet_face.ver & 3][tet_face.tet].emplace_back(i);
            const TriFace& nei = mesh.tets[tet_face.tet].nei[tet_face.ver & 3];
            tet_map[nei.ver & 3][nei.tet].emplace_back(i);
            continue;
        }
    }
}
