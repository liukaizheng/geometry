#include <graphcut/graphcut.h>
#include <polygonalization/bsp_complex.h>
#include <utils/disjoint-set.h>

#include <Eigen/Dense>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <queue>
#include <tuple>
#include <unordered_map>

inline uint32_t remove_ghost_tets(const TetMesh& mesh, std::vector<uint32_t>& new_order) {
    const auto& tets = mesh.tets;
    uint32_t idx = 0;
    new_order.resize(tets.size(), TriFace::INVALID);
    for (uint32_t i = 0; i < tets.size(); i++) {
        if (!mesh.is_hull_tet(i)) {
            new_order[i] = idx++;
        }
    }
    return idx;
}

inline bool tet_face_is_new(const uint32_t tet_ind, const uint32_t adj_tet_ind, const uint32_t adj_cell_ind) {
    return tet_ind < adj_tet_ind || adj_cell_ind == TriFace::INVALID;
}

inline uint32_t find_tet_edge(
    uint32_t ev1, uint32_t ev2, std::unordered_map<uint32_t, uint32_t>* edge_map, std::vector<BSPEdge>& edges
) {
    if (ev1 > ev2) {
        std::swap(ev1, ev2);
    }
    auto& map = edge_map[ev1];
    const auto it = map.find(ev2);
    if (it != map.end()) {
        return it->second;
    } else {
        const uint32_t eid = static_cast<uint32_t>(edges.size());
        map.emplace(ev2, eid);
        edges.emplace_back(ev1, ev2, eid);
        return eid;
    }
}

inline void fill_face_color(
    const uint32_t tid, const uint32_t fid, const std::vector<std::vector<uint32_t>>& constraints, BSPComplex* complex
) {
    if (constraints[tid].empty()) {
        complex->faces[fid].color = FaceColor::WHITE;
    } else {
        BSPFace& face = complex->faces[fid];
        for (const uint32_t cid : constraints[tid]) {
            if (cid < complex->constraints->n_triangles) {
                face.coplanar_constraints.emplace_back(cid);
            }
        }
        if (face.coplanar_constraints.empty()) {
            face.color = FaceColor::WHITE;
        } else {
            face.color = FaceColor::GRAY;
        }
    }
}

BSPComplex::BSPComplex(
    const TetMesh& mesh, const Constraints* c,
    const std::vector<std::unordered_map<uint32_t, uint32_t>>& ori_edge_parents,
    std::array<std::vector<std::vector<uint32_t>>, 5>&& tet_maps
)
    : constraints(c) {
    vertices.reserve(mesh.n_points);
    for (uint32_t i = 0; i < mesh.n_points; i++) {
        const double* p = mesh.point(i);
        vertices.emplace_back(new ExplicitPoint3D(p[0], p[1], p[2]));
    }

    verts_oris.resize(vertices.size(), 2);
    std::vector<uint32_t> new_order;
    const uint32_t n_cells = remove_ghost_tets(mesh, new_order);
    cells.resize(n_cells);
    edges.reserve(n_cells + mesh.n_points);
    faces.reserve(cells.size() << 1);

    std::vector<std::unordered_map<uint32_t, uint32_t>> edge_map(mesh.n_points - 1);
    auto push_edge = [&edge_map, &ori_edge_parents, this](uint32_t pa, uint32_t pb) {
        if (pa > pb) {
            std::swap(pa, pb);
        }
        if (edge_map[pa].find(pb) != edge_map[pa].end()) {
            return;
        }
        const auto& map = ori_edge_parents[pa];
        if (map.find(pb) != map.end()) {
            const uint32_t eid = static_cast<uint32_t>(edges.size());
            edges.emplace_back(pa, pb, eid);
            edge_map[pa].emplace(pb, eid);
        }
    };
    for (uint32_t i = 0; i < mesh.tets.size(); i++) {
        const uint32_t cell_ind = new_order[i];
        if (cell_ind == TriFace::INVALID) {
            continue;
        }
        const Tet& tet = mesh.tets[i];
        for (uint32_t j = 0; j < 3; j++) {
            const uint32_t va = tet.data[j];
            for (uint32_t k = j + 1; k < 4; k++) {
                push_edge(va, tet.data[k]);
            }
        }
    }
    n_ori_edges = static_cast<uint32_t>(edges.size());

    std::vector<uint32_t> tet_edges;
    tet_edges.reserve(6);
    for (uint32_t i = 0; i < mesh.tets.size(); i++) {
        const uint32_t cell_ind = new_order[i];
        if (cell_ind == TriFace::INVALID) {
            continue;
        }
        cells[cell_ind].constraints = tet_maps[4][i];
        constexpr std::array<uint32_t, 12> edgetbl{{3, 5, 4, 1, 5, 2, 2, 4, 0, 0, 3, 1}};
        const Tet& tet = mesh.tets[i];
        tet_edges.clear();
        for (uint32_t j = 0; j < 3; j++) {
            const uint32_t va = tet.data[j];
            for (uint32_t k = j + 1; k < 4; k++) {
                tet_edges.emplace_back(find_tet_edge(va, tet.data[k], edge_map.data(), edges));
            }
        }
        for (uint32_t j = 0; j < 4; j++) {
            const auto& nei = tet.nei[j];
            const uint32_t adj_cell_ind = new_order[nei.tet];
            if (tet_face_is_new(i, nei.tet, adj_cell_ind)) {
                const uint32_t verts[]{mesh.org(nei), mesh.dest(nei), mesh.apex(nei)};
                const uint32_t fid = static_cast<uint32_t>(faces.size());
                faces.emplace_back(verts, fid, cell_ind, adj_cell_ind);
                cells[cell_ind].faces.emplace_back(fid);
                if (adj_cell_ind != TriFace::INVALID) {
                    cells[adj_cell_ind].faces.emplace_back(fid);
                }
                const uint32_t* edge_indices = &edgetbl[j * 3];
                for (uint32_t k = 0; k < 3; k++) {
                    const uint32_t eid = tet_edges[edge_indices[k]];
                    faces.back().edges.emplace_back(eid);
                    edges[eid].face = fid;
                }
                fill_face_color(i, fid, tet_maps[j], this);
            }
        }
    }
    vert_visit.resize(vertices.size(), 0);
    edge_visit.resize(edges.size(), 0);
}

inline bool is_point_built_from_plane(
    const GenericPoint3D* p, const ExplicitPoint3D* p0, const ExplicitPoint3D* p1, const ExplicitPoint3D* p2
) {
    if (p->is_explicit()) {
        return p == p0 || p == p1 || p == p2;
    } else if (p->is_lpi()) {
        const ImplicitPointLPI& lpi = p->to_lpi();
        return (
            (((p0 == &lpi.p()) && (p1 == &lpi.q())) || ((p1 == &lpi.p()) && (p0 == &lpi.q()))) ||
            (((p1 == &lpi.p()) && (p2 == &lpi.q())) || ((p2 == &lpi.p()) && (p1 == &lpi.q()))) ||
            (((p2 == &lpi.p()) && (p0 == &lpi.q())) || ((p0 == &lpi.p()) && (p2 == &lpi.q()))) ||
            ((p0 == &lpi.r()) && (p1 == &lpi.s()) && (p2 == &lpi.t()))
        );
    } else {
        const ImplicitPointTPI& tpi = p->to_tpi();
        return (
            ((p0 == &tpi.v1()) && (p1 == &tpi.v2()) && (p2 == &tpi.v3())) ||
            ((p0 == &tpi.w1()) && (p1 == &tpi.w2()) && (p2 == &tpi.w3())) ||
            ((p0 == &tpi.u1()) && (p1 == &tpi.u2()) && (p2 == &tpi.u3()))
        );
    }
}

inline void
verts_orient_wrt_plane(BSPComplex* complex, const uint32_t* c, const uint32_t* verts, const uint32_t n_verts) {
    const auto& p0 = complex->vertices[c[0]]->to_explicit();
    const auto& p1 = complex->vertices[c[1]]->to_explicit();
    const auto& p2 = complex->vertices[c[2]]->to_explicit();
    for (uint32_t i = 0; i < n_verts; i++) {
        const uint32_t vid = verts[i];
        if (complex->verts_oris[vid] != 2) {
            continue;
        }
        const auto* p = complex->vertices[vid];
        if (is_point_built_from_plane(p, &p0, &p1, &p2)) {
            complex->verts_oris[vid] = 0;
        } else {
            complex->verts_oris[vid] = GenericPoint3D::orient3d(*p, p0, p2, p1);
        }
    }
}

inline void find_coplanar_constraints(
    BSPComplex* complex, std::vector<uint32_t>& cell_constraints, const uint32_t* c, std::vector<uint32_t>& coplanar_c
) {
    const auto it = std::remove_if(cell_constraints.begin(), cell_constraints.end(), [complex, c](const uint32_t cid) {
        if (complex->constraints->is_virtual(cid)) {
            return false;
        }
        const uint32_t* tri = &complex->constraints->triangles[cid * 3];
        verts_orient_wrt_plane(complex, c, tri, 3);
        return complex->verts_oris[tri[0]] == 0 && complex->verts_oris[tri[1]] == 0 && complex->verts_oris[tri[2]] == 0;
    });
    std::move(it, cell_constraints.end(), std::back_inserter(coplanar_c));
    cell_constraints.resize(static_cast<uint32_t>(std::distance(cell_constraints.begin(), it)));
}

inline void find_cell_verts_and_edges(
    BSPComplex* complex, const BSPCell& cell, std::vector<uint32_t>& cell_verts, std::vector<uint32_t>& cell_edges
) {
    std::vector<uint32_t>& vert_visit = complex->vert_visit;
    std::vector<int>& edge_visit = complex->edge_visit;
    for (uint32_t i = 0; i < cell.faces.size(); i++) {
        const BSPFace& face = complex->faces[cell.faces[i]];
        for (const uint32_t eid : face.edges) {
            const BSPEdge& edge = complex->edges[eid];
            if (edge_visit[eid] == 0) {
                edge_visit[eid] = 1;
                cell_edges.emplace_back(eid);
                for (const uint32_t ev : edge.vertices) {
                    if (vert_visit[ev] == 0) {
                        vert_visit[ev] = 1;
                        cell_verts.emplace_back(ev);
                    }
                }
            }
        }
    }
    for (const uint32_t v : cell_verts) {
        vert_visit[v] = 0;
    }
    for (const uint32_t e : cell_edges) {
        edge_visit[e] = 0;
    }
}

inline void count_vert_orient(
    const std::vector<uint32_t>& verts, const int* orient, uint32_t& n_over, uint32_t& n_under, uint32_t& n_on
) {
    for (const uint32_t v : verts) {
        if (orient[v] == 0) {
            n_on += 1;
        } else if (orient[v] == 1) {
            n_over += 1;
        } else if (orient[v] == -1) {
            n_under += 1;
        }
    }
}

inline bool constraint_inner_intersects_edge(const BSPEdge& edge, const int* oris) {
    return (oris[edge.vertices[0]] > 0 && oris[edge.vertices[1]] < 0) ||
           (oris[edge.vertices[0]] < 0 && oris[edge.vertices[1]] > 0);
}

inline uint32_t oppo_edge_face(
    const std::vector<BSPFace>& faces, const uint32_t eid, const uint32_t fid, const std::vector<uint32_t>& cell_faces
) {
    for (const uint32_t f : cell_faces) {
        if (f == fid) {
            continue;
        }
        const auto& face_edges = faces[f].edges;
        const auto it = std::find(face_edges.begin(), face_edges.end(), eid);
        if (it != face_edges.end()) {
            return f;
        }
    }
    return TriFace::INVALID;
}

inline uint32_t oppo_cell(const uint32_t cid, const uint32_t* cells) { return cells[0] == cid ? cells[1] : cells[0]; }

inline uint32_t add_lpi_vert(BSPComplex* complex, const BSPEdge& edge, const uint32_t cid) {
    const uint32_t* tri = &complex->constraints->triangles[cid * 3];
    auto& vertices = complex->vertices;
    const uint32_t vid = static_cast<uint32_t>(vertices.size());
    complex->vertices.emplace_back(new ImplicitPointLPI(
        vertices[edge.mesh_vertices[0]]->to_explicit(), vertices[edge.mesh_vertices[1]]->to_explicit(),
        vertices[tri[0]]->to_explicit(), vertices[tri[1]]->to_explicit(), vertices[tri[2]]->to_explicit()
    ));
    complex->verts_oris.emplace_back(2);
    complex->vert_visit.emplace_back(0);
    return vid;
}

inline bool two_equal_vertices(const uint32_t* v1, const uint32_t* v2, uint32_t* comm) {
    uint32_t i = 0;
    if (v1[0] == v2[0] || v1[0] == v2[1] || v1[0] == v2[2]) {
        comm[i++] = v1[0];
    }
    if (v1[1] == v2[0] || v1[1] == v2[1] || v1[1] == v2[2]) {
        comm[i++] = v1[1];
    }
    if (v1[2] == v2[0] || v1[2] == v2[1] || v1[2] == v2[2]) {
        comm[i++] = v1[2];
    }
    return i >= 2;
}

inline uint32_t add_tpi_vert(BSPComplex* complex, const BSPEdge& edge, const uint32_t cid) {
    const uint32_t* tri = &complex->constraints->triangles[cid * 3];
    auto& vertices = complex->vertices;
    const uint32_t vid = static_cast<uint32_t>(vertices.size());
    uint32_t comm[3];
    // clang-format off
    if (two_equal_vertices(tri, edge.mesh_vertices, comm)) {
        complex->vertices.emplace_back(new ImplicitPointLPI(
                vertices[comm[0]]->to_explicit(),
                vertices[comm[1]]->to_explicit(),
                vertices[edge.mesh_vertices[3]]->to_explicit(),
                vertices[edge.mesh_vertices[4]]->to_explicit(),
                vertices[edge.mesh_vertices[5]]->to_explicit()
        ));
    } else if (two_equal_vertices(tri, &edge.mesh_vertices[3], comm)) {
        complex->vertices.emplace_back(new ImplicitPointLPI(
                vertices[comm[0]]->to_explicit(),
                vertices[comm[1]]->to_explicit(),
                vertices[edge.mesh_vertices[0]]->to_explicit(),
                vertices[edge.mesh_vertices[1]]->to_explicit(),
                vertices[edge.mesh_vertices[2]]->to_explicit()
        ));
    } else if (two_equal_vertices(edge.mesh_vertices, &edge.mesh_vertices[3], comm)) {
        complex->vertices.emplace_back(new ImplicitPointLPI(
                vertices[comm[0]]->to_explicit(),
                vertices[comm[1]]->to_explicit(),
                vertices[tri[0]]->to_explicit(),
                vertices[tri[1]]->to_explicit(),
                vertices[tri[2]]->to_explicit()
        ));
    } else {
        complex->vertices.emplace_back(new ImplicitPointTPI(
                vertices[edge.mesh_vertices[0]]->to_explicit(),
                vertices[edge.mesh_vertices[1]]->to_explicit(),
                vertices[edge.mesh_vertices[2]]->to_explicit(),
                vertices[edge.mesh_vertices[3]]->to_explicit(),
                vertices[edge.mesh_vertices[4]]->to_explicit(),
                vertices[edge.mesh_vertices[5]]->to_explicit(),
                vertices[tri[0]]->to_explicit(),
                vertices[tri[1]]->to_explicit(),
                vertices[tri[2]]->to_explicit()
        ));
    }
    // clang-format on
    complex->verts_oris.emplace_back(2);
    complex->vert_visit.emplace_back(0);
    return vid;
}

inline void incident(const BSPComplex* complex, const uint32_t eid, std::vector<uint32_t>& ef) {
    const auto& edge = complex->edges[eid];
    uint32_t f = edge.face;
    uint32_t c = complex->faces[f].cells[0];
    ef.emplace_back(f);
    while (true) {
        f = oppo_edge_face(complex->faces, eid, f, complex->cells[c].faces);
        if (f == edge.face) {
            return;
        } else {
            ef.emplace_back(f);
        }
        c = oppo_cell(c, complex->faces[f].cells);
        if (c == TriFace::INVALID) {
            break;
        }
    }
    f = edge.face;
    if ((c = complex->faces[f].cells[1]) == TriFace::INVALID) {
        return;
    }
    while (true) {
        f = oppo_edge_face(complex->faces, eid, f, complex->cells[c].faces);
        ef.emplace_back(f);
        c = oppo_cell(c, complex->faces[f].cells);
        if (c == TriFace::INVALID) {
            break;
        }
    }
}

inline bool consecutive_edges(const uint32_t* ev1, const uint32_t* ev2) {
    return ev1[0] == ev2[0] || ev1[0] == ev2[1] || ev1[1] == ev2[0] || ev1[1] == ev2[1];
}

inline void add_edge_to_face(BSPComplex* complex, const uint32_t eid, std::vector<uint32_t>& face_edges) {
    BSPEdge& new_edge = complex->edges[eid];
    const uint32_t n_edges = static_cast<uint32_t>(face_edges.size());
    for (uint32_t i = 0; i < n_edges; i++) {
        if (consecutive_edges(new_edge.vertices.data(), complex->edges[face_edges[i]].vertices.data())) {
            const uint32_t pos = (i + 1) % n_edges;
            if (consecutive_edges(new_edge.vertices.data(), complex->edges[face_edges[pos]].vertices.data())) {
                face_edges.emplace(face_edges.begin() + pos, eid);
            } else {
                face_edges.emplace(face_edges.begin() + i, eid);
            }
            return;
        }
    }
}

inline void split_edge(BSPComplex* complex, const uint32_t eid, const uint32_t constr_id) {
    BSPEdge& edge = complex->edges[eid];
    std::vector<uint32_t> ef;
    incident(complex, eid, ef);
    const uint32_t new_vid = edge.mesh_vertices[2] == TriFace::INVALID ? add_lpi_vert(complex, edge, constr_id)
                                                                       : add_tpi_vert(complex, edge, constr_id);
    const uint32_t new_eid = static_cast<uint32_t>(complex->edges.size());
    complex->edges.emplace_back(edge.split(new_vid));
    complex->edge_visit.emplace_back(0);
    for (const uint32_t f : ef) {
        add_edge_to_face(complex, new_eid, complex->faces[f].edges);
    }
}

inline void face_vertices(const BSPComplex* complex, const BSPFace& face, std::vector<uint32_t>& vertices) {
    vertices.reserve(face.edges.size());
    const uint32_t* first = complex->edges[face.edges[0]].vertices.data();
    const uint32_t* last = complex->edges[face.edges.back()].vertices.data();
    if (first[0] == last[0] || first[0] == last[1]) {
        vertices.emplace_back(first[0]);
        vertices.emplace_back(first[1]);
    } else {
        vertices.emplace_back(first[1]);
        vertices.emplace_back(first[0]);
    }
    for (uint32_t i = 2; i < face.edges.size(); i++) {
        const BSPEdge& edge = complex->edges[face.edges[i - 1]];
        if (edge.vertices[0] == vertices.back()) {
            vertices.emplace_back(edge.vertices[1]);
        } else {
            vertices.emplace_back(edge.vertices[0]);
        }
    }
}

inline bool constraint_inner_intersects_face(const std::vector<uint32_t>& face_verts, int* vert_orient) {
    uint32_t n_over = 0, n_under = 0, n_on = 0;
    count_vert_orient(face_verts, vert_orient, n_over, n_under, n_on);
    return n_over > 0 && n_under > 0;
}

inline void edges_partition(BSPComplex* complex, const uint32_t fid, const uint32_t new_fid) {
    BSPFace& face = complex->faces[fid];
    uint32_t idx = 0;
    const uint32_t n_edges = static_cast<uint32_t>(face.edges.size());
    const auto& vert_oris = complex->verts_oris;
    // search the first edge having first vertex == 0 and second vertex < 0
    for (; idx < face.edges.size(); idx++) {
        const uint32_t eid = face.edges[idx];
        const uint32_t next_eid = face.edges[(idx + 1) % n_edges];
        const BSPEdge& edge = complex->edges[eid];
        const BSPEdge& next_edge = complex->edges[next_eid];
        const uint32_t comm_vert =
            edge.vertices[0] == next_edge.vertices[0] || edge.vertices[0] == next_edge.vertices[1] ? 0 : 1;
        if (vert_oris[edge.vertices[comm_vert]] < 0 && vert_oris[edge.vertices[!comm_vert]] == 0) {
            break;
        }
    }
    std::rotate(face.edges.begin(), face.edges.begin() + idx, face.edges.end());
    for (idx = 1; idx < face.edges.size(); idx++) {
        const BSPEdge& edge = complex->edges[face.edges[idx]];
        if (vert_oris[edge.vertices[0]] == 0 || vert_oris[edge.vertices[1]] == 0) {
            break;
        }
    }
    idx += 1;
    complex->faces[new_fid].edges.assign(face.edges.begin() + idx, face.edges.end());
    face.edges.resize(idx);
    for (const uint32_t eid : complex->faces[new_fid].edges) {
        complex->edges[eid].face = new_fid;
    }
}

inline void add_common_edge(
    BSPComplex* complex, const uint32_t fid, const uint32_t new_fid, const uint32_t constr_id, const uint32_t* endpts
) {
    BSPFace& face = complex->faces[fid];
    BSPFace& new_face = complex->faces[new_fid];
    const uint32_t new_eid = static_cast<uint32_t>(complex->edges.size());
    complex->edges.emplace_back(
        endpts[0], endpts[1], new_eid, face.mesh_vertices, &complex->constraints->triangles[constr_id * 3]
    );
    complex->edge_visit.emplace_back(0);
    complex->edges.back().face = fid;
    face.edges.emplace_back(new_eid);
    new_face.edges.emplace_back(new_eid);
}

inline void
split_face(BSPComplex* complex, const uint32_t fid, const uint32_t constr_id, const std::vector<uint32_t>& face_verts) {
    const BSPFace& face = complex->faces[fid];
    const uint32_t new_fid = static_cast<uint32_t>(complex->faces.size());
    complex->faces.emplace_back(face.mesh_vertices, face.parent, face.cells[0], face.cells[1]);
    complex->faces.back().color = face.color;
    complex->faces.back().coplanar_constraints = face.coplanar_constraints;

    const auto& face_cells = complex->faces[fid].cells;
    complex->cells[face_cells[0]].faces.emplace_back(new_fid);
    if (face_cells[1] != TriFace::INVALID) {
        complex->cells[face_cells[1]].faces.emplace_back(new_fid);
    }

    uint32_t zero_verts[2];
    uint32_t pos = 0;
    for (const uint32_t vid : face_verts) {
        if (complex->verts_oris[vid] == 0) {
            zero_verts[pos++] = vid;
        }
    }

    edges_partition(complex, fid, new_fid);
    add_common_edge(complex, fid, new_fid, constr_id, zero_verts);
}

inline void move_face(BSPComplex* complex, const uint32_t cell_fid, const uint32_t cid, const uint32_t new_cid) {
    const uint32_t fid = complex->cells[cid].faces[cell_fid];
    BSPFace& face = complex->faces[fid];
    if (face.cells[0] == cid) {
        face.cells[0] = new_cid;
    } else {
        face.cells[1] = new_cid;
    }
    complex->cells[new_cid].faces.emplace_back(fid);
    complex->cells[cid].remove_face(cell_fid);
}

inline void faces_partition(BSPComplex* complex, const uint32_t cid, const uint32_t new_cid) {
    BSPCell& cell = complex->cells[cid];
    uint32_t n_faces = static_cast<uint32_t>(cell.faces.size());
    uint32_t idx = 0;
    while (idx < n_faces) {
        const uint32_t fid = cell.faces[idx];
        const BSPFace& face = complex->faces[fid];
        std::vector<uint32_t> face_verts;
        face_vertices(complex, face, face_verts);
        bool found = false;
        for (const uint32_t vid : face_verts) {
            if (complex->verts_oris[vid] > 0) {
                move_face(complex, idx, cid, new_cid);
                found = true;
                n_faces -= 1;
                break;
            }
        }
        if (!found) {
            idx += 1;
        }
    }
}

inline void add_edges_to_face(BSPComplex* complex, BSPFace& face, const std::vector<uint32_t>& edges) {
    std::vector<uint32_t> face_verts(edges.size());
    auto& vert_visit = complex->vert_visit;
    uint32_t idx = 0;
    for (const uint32_t eid : edges) {
        const BSPEdge& edge = complex->edges[eid];
        for (const uint32_t vid : edge.vertices) {
            if (vert_visit[vid] == 0) {
                face_verts[idx++] = vid;
                vert_visit[vid] = 1;
            }
        }
    }

    for (idx = 0; idx < face_verts.size(); idx++) {
        vert_visit[face_verts[idx]] = idx;
    }

    std::vector<uint32_t> vert_edges(face_verts.size() << 1, TriFace::INVALID);
    for (const uint32_t eid : edges) {
        for (const uint32_t vid : complex->edges[eid].vertices) {
            uint32_t* v_edges = &vert_edges[vert_visit[vid] << 1];
            if (v_edges[0] == TriFace::INVALID) {
                v_edges[0] = eid;
            } else {
                v_edges[1] = eid;
            }
        }
    }

    face.edges.resize(edges.size());

    idx = 0;
    uint32_t nv = face_verts[0];
    uint32_t e = vert_edges[vert_visit[face_verts[0]] << 1];
    do {
        face.edges[idx++] = e;
        const uint32_t* v_edges = &vert_edges[vert_visit[nv] << 1];
        e = v_edges[0] == e ? v_edges[1] : v_edges[0];
        const uint32_t* evs = complex->edges[e].vertices.data();
        nv = evs[0] == nv ? evs[1] : evs[0];
    } while (idx < edges.size());

    for (const uint32_t vid : face_verts) {
        vert_visit[vid] = 0;
    }
}

inline bool edges_share_common_plane(const BSPEdge& a, const BSPEdge& b) {
    const uint32_t* mva = a.mesh_vertices;
    const uint32_t* mvb = b.mesh_vertices;
    return (
        mva[0] == mvb[0] && mva[1] == mvb[1] && mva[2] == mvb[2] && mva[3] == mvb[3] && mva[4] == mvb[4] &&
        mva[5] == mvb[5]
    );
}

inline void fix_common_face_orientation(BSPComplex* complex, const uint32_t fid) {
    BSPFace& face = complex->faces[fid];
    const uint32_t* mv = face.mesh_vertices;
    double mvc[9]; // mesh vertex coordinates
    complex->vertices[mv[0]]->to_double_approx(mvc);
    complex->vertices[mv[1]]->to_double_approx(mvc + 3);
    complex->vertices[mv[2]]->to_double_approx(mvc + 6);
    const int xyz = GenericPoint3D::max_component_at_triangle_normal(mvc, mvc + 3, mvc + 6);
    double d_base_ori = 0.0;
    if (xyz == 2) {
        d_base_ori = orient2d(mvc, mvc + 3, mvc + 6);
    } else if (xyz == 0) {
        d_base_ori = orient2d(mvc + 1, mvc + 4, mvc + 7);
    } else {
        double p0[2]{mvc[2], mvc[0]};
        double p1[2]{mvc[5], mvc[3]};
        double p2[2]{mvc[8], mvc[6]};
        d_base_ori = orient2d(p0, p1, p2);
    }
    int base_ori = (d_base_ori > 0.0) - (d_base_ori < 0.0);

    const BSPEdge& edge0 = complex->edges[face.edges[0]];
    const uint32_t* first = edge0.vertices.data();
    const uint32_t* last = complex->edges[face.edges.back()].vertices.data();
    uint32_t v0, v1;
    if (first[0] == last[0] || first[0] == last[1]) {
        v0 = first[0];
        v1 = first[1];
    } else {
        v0 = first[1];
        v1 = first[0];
    }
    for (uint32_t i = 1; i < face.edges.size(); i++) {
        const BSPEdge& edge = complex->edges[face.edges[i]];
        if (edges_share_common_plane(edge0, edge)) {
            continue;
        }
        const uint32_t v2 = edge.vertices[0] == v1 ? edge.vertices[1] : edge.vertices[0];
        int ori = base_ori *
                  GenericPoint3D::orient2d(*complex->vertices[v0], *complex->vertices[v1], *complex->vertices[v2], xyz);
        if (ori < 0) {
            return;
        } else {
            std::swap(face.cells[0], face.cells[1]);
        }
    }
}

inline void add_common_face(
    BSPComplex* complex, const uint32_t constr_id, const uint32_t cid, const uint32_t new_cid,
    const std::vector<uint32_t>& cell_edges
) {
    const uint32_t new_fid = static_cast<uint32_t>(complex->faces.size());
    complex->faces.emplace_back(&complex->constraints->triangles[constr_id * 3], new_fid, cid, new_cid);
    complex->faces.back().color = complex->constraints->is_virtual(constr_id) ? FaceColor::WHITE : FaceColor::GRAY;
    std::vector<uint32_t> comm_edges;
    const auto& vert_oris = complex->verts_oris;
    for (const uint32_t eid : cell_edges) {
        const BSPEdge& edge = complex->edges[eid];
        if (vert_oris[edge.vertices[0]] == 0 && vert_oris[edge.vertices[1]] == 0) {
            comm_edges.emplace_back(eid);
        }
    }

    add_edges_to_face(complex, complex->faces[new_fid], comm_edges);
    for (const uint32_t eid : comm_edges) {
        complex->edges[eid].face = new_fid;
    }
    complex->cells[cid].faces.emplace_back(new_fid);
    complex->cells[new_cid].faces.emplace_back(new_fid);
    fix_common_face_orientation(complex, new_fid);
}

inline void remove_ele(std::vector<uint32_t>& vec, const uint32_t ind) {
    if (ind + 1 != vec.size()) {
        vec[ind] = vec.back();
    }
    vec.pop_back();
}

inline void
constraints_partition(BSPComplex* complex, const uint32_t constr_id, const uint32_t down_cid, const uint32_t up_cid) {
    BSPCell& down_cell = complex->cells[down_cid];
    BSPCell& up_cell = complex->cells[up_cid];
    const uint32_t* constraint = &complex->constraints->triangles[constr_id * 3];
    uint32_t n_constr = static_cast<uint32_t>(down_cell.constraints.size());
    auto& vert_oris = complex->verts_oris;
    for (uint32_t i = 0; i < n_constr;) {
        const uint32_t* triangle = &complex->constraints->triangles[down_cell.constraints[i] * 3];
        verts_orient_wrt_plane(complex, constraint, triangle, 3);
        uint32_t n_over = 0, n_under = 0;
        for (uint32_t j = 0; j < 3; j++) {
            if (vert_oris[triangle[j]] == 1) {
                n_over += 1;
            } else if (vert_oris[triangle[j]] == -1) {
                n_under += 1;
            }
        }
        if (n_over == 0 && n_under == 0) { // this constraint is virtual
            vert_oris[triangle[0]] = 2;
            vert_oris[triangle[1]] = 2;
            vert_oris[triangle[2]] = 2;
            remove_ele(down_cell.constraints, i);
            n_constr -= 1;
            continue;
        }
        if (n_over > 0) {
            up_cell.constraints.emplace_back(down_cell.constraints[i]);
        }
        if (n_under == 0) {
            remove_ele(down_cell.constraints, i);
            n_constr -= 1;
            continue;
        }
        i += 1;
    }
}

void BSPComplex::split_cell(const uint32_t cid) {
    /*{
        double nv[3];
        const auto vv = vertices.back();
        vv->to_double(nv);
        if (nv[2] > 0) {
            const int b = 2;
        }
        BSPCell& _c = cells[172];
        std::vector<uint32_t> v, e;
        find_cell_verts_and_edges(this, _c, v, e);
        std::vector<double> p(v.size() * 3);
        for (uint32_t i = 0; i < v.size(); i++) {
            vertices[v[i]]->to_double(&p[i * 3]);
        }
        std::vector<std::vector<uint32_t>> evs;
        for (const uint32_t eid : e) {
            evs.emplace_back(std::vector{edges[eid].vertices[0], edges[eid].vertices[1]});
        }
        std::vector<std::vector<uint32_t>> fvs;
        for (const auto fid : _c.faces) {
            std::vector<uint32_t> fv;
            face_vertices(this, faces[fid], fv);
            fvs.emplace_back(fv);
        }
        const auto* p1 = vertices[1];
        const auto* p2 = vertices[15];
        const auto* p3 = vertices[10];
        std::vector<int> oo;
        for (const uint32_t vid : v) {
            oo.emplace_back(GenericPoint3D::orient3d(*vertices[vid], *p1, *p2, *p3));
        }
        const int a = 2;
    }*/
    BSPCell& cell = cells[cid];
    const uint32_t constr_id = cell.constraints.back();
    cell.constraints.pop_back();
    const uint32_t* constraint = &constraints->triangles[constr_id * 3];
    std::vector<uint32_t> coplanar_c;
    find_coplanar_constraints(this, cell.constraints, constraint, coplanar_c);
    if (!constraints->is_virtual(constr_id)) {
        coplanar_c.emplace_back(constr_id);
    }

    std::vector<uint32_t> cell_verts;
    std::vector<uint32_t> cell_edges;
    find_cell_verts_and_edges(this, cell, cell_verts, cell_edges);

    verts_orient_wrt_plane(this, constraint, cell_verts.data(), static_cast<uint32_t>(cell_verts.size()));

    uint32_t n_over = 0, n_under = 0, n_on = 0;
    count_vert_orient(cell_verts, verts_oris.data(), n_over, n_under, n_on);
    if (n_over == 0 || n_under == 0) {
        std::fill(verts_oris.begin(), verts_oris.end(), 2);
        return;
    }

    for (uint32_t i = 0, ne = static_cast<uint32_t>(cell_edges.size()); i < ne; i++) {
        const uint32_t eid = cell_edges[i];
        const BSPEdge& edge = edges[eid];
        if (constraint_inner_intersects_edge(edge, verts_oris.data())) {
            split_edge(this, eid, constr_id);
            const uint32_t new_vert = static_cast<uint32_t>(vertices.size() - 1);
            cell_verts.emplace_back(new_vert);
            cell_edges.emplace_back(static_cast<uint32_t>(edges.size() - 1));
            verts_oris[new_vert] = 0;
        }
    }

    const auto n_faces = cell.faces.size();
    for (uint32_t i = 0; i < n_faces; i++) {
        const uint32_t fid = cell.faces[i];
        BSPFace& face = faces[fid];
        std::vector<uint32_t> face_verts;
        face_vertices(this, face, face_verts);
        if (constraint_inner_intersects_face(face_verts, verts_oris.data())) {
            split_face(this, fid, constr_id, face_verts);
            cell_edges.emplace_back(edges.size() - 1);
        }
    }

    const uint32_t new_cid = static_cast<uint32_t>(cells.size());
    cells.emplace_back();
    faces_partition(this, cid, new_cid);
    add_common_face(this, constr_id, cid, new_cid, cell_edges);
    faces.back().coplanar_constraints = std::move(coplanar_c);
    constraints_partition(this, constr_id, cid, new_cid);

    // reset vertex orientation wrt constraint
    std::fill(verts_oris.begin(), verts_oris.end(), 2);
}

inline int face_normal_dominant_component(const BSPComplex* complex, const BSPFace& face) {
    const uint32_t* mv = face.mesh_vertices;
    double mvc[9]; // mesh vertex coordinates
    complex->vertices[mv[0]]->to_double_approx(mvc);
    complex->vertices[mv[1]]->to_double_approx(mvc + 3);
    complex->vertices[mv[2]]->to_double_approx(mvc + 6);
    return GenericPoint3D::max_component_at_triangle_normal(mvc, mvc + 3, mvc + 6);
}

inline void face_baricenter_approx(const BSPComplex* complex, const BSPFace& face, double* bar) {
    bar[0] = 0.0;
    bar[1] = 0.0;
    bar[2] = 0.0;

    const uint32_t* last = complex->edges[face.edges.back()].vertices.data();
    const uint32_t* first = complex->edges[face.edges[0]].vertices.data();
    uint32_t vid = first[0] == last[0] || first[0] == last[1] ? first[0] : first[1];
    double p[3];
    for (const uint32_t eid : face.edges) {
        const uint32_t* e_verts = complex->edges[eid].vertices.data();
        vid = e_verts[0] == vid ? e_verts[1] : e_verts[0];
        complex->vertices[vid]->to_double_approx(p);
        bar[0] += p[0];
        bar[1] += p[1];
        bar[2] += p[2];
    }
    bar[0] /= static_cast<double>(face.edges.size());
    bar[1] /= static_cast<double>(face.edges.size());
    bar[2] /= static_cast<double>(face.edges.size());
}

inline bool is_baricenter_on_face(
    const BSPComplex* complex, const BSPFace& face, const ExplicitPoint3D& face_center, const int xyz
) {
    const uint32_t* last = complex->edges[face.edges.back()].vertices.data();
    const uint32_t* first = complex->edges[face.edges[0]].vertices.data();
    uint32_t vid = first[0] == last[0] || first[0] == last[1] ? first[0] : first[1];
    int base_ori = 0;
    for (const uint32_t eid : face.edges) {
        const uint32_t* e_verts = complex->edges[eid].vertices.data();
        const uint32_t pvid = vid;
        vid = e_verts[0] == vid ? e_verts[1] : e_verts[0];
        const int ori = GenericPoint3D::orient2d(face_center, *complex->vertices[pvid], *complex->vertices[vid], xyz);
        if (ori == 0) {
            return false;
        }
        if (ori != base_ori) {
            if (base_ori == 0) {
                base_ori = ori;
            } else {
                return false;
            }
        }
    }
    return true;
}

// return 1 if point is on the boundary of the triangle, 2 if it is in the interior, 0 otherwise
inline int localized_point_on_triangle(
    const GenericPoint3D& p, const GenericPoint3D& a, const GenericPoint3D& b, const GenericPoint3D& c, int xyz
) {
    int o1, o2, o3;
    if (xyz == 2) {
        o1 = GenericPoint3D::orient_xy(p, a, b);
        o2 = GenericPoint3D::orient_xy(p, b, c);
        o3 = GenericPoint3D::orient_xy(p, c, a);
    } else if (xyz == 0) {
        o1 = GenericPoint3D::orient_yz(p, a, b);
        o2 = GenericPoint3D::orient_yz(p, b, c);
        o3 = GenericPoint3D::orient_yz(p, c, a);
    } else {
        o1 = GenericPoint3D::orient_zx(p, a, b);
        o2 = GenericPoint3D::orient_zx(p, b, c);
        o3 = GenericPoint3D::orient_zx(p, c, a);
    }
    return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0)) +
           ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
}

inline bool coplanar_constraint_inner_intersects_face(
    const BSPComplex* complex, const std::vector<uint32_t>& f_edges, const uint32_t* tri, const int xyz
) {
    int mask = 0;
    const uint32_t* last = complex->edges[f_edges.back()].vertices.data();
    const uint32_t* first = complex->edges[f_edges[0]].vertices.data();
    const uint32_t vid_0 = first[0] == last[0] || first[0] == last[1] ? first[0] : first[1];
    for (int i = 0; i < 3; i++) {
        uint32_t vid = vid_0;
        for (const uint32_t eid : f_edges) {
            const uint32_t* e_verts = complex->edges[eid].vertices.data();
            vid = e_verts[0] == vid ? e_verts[1] : e_verts[0];
            if (tri[i] == vid) {
                mask |= (1 << i);
                break;
            }
        }
    }
    if (mask == 7) return true;
    const auto& vertices = complex->vertices;
    for (int i = 0; i < 3; i++) {
        if ((mask & (1 << i)) != 0) {
            continue;
        }
        for (uint32_t e = 0; e < f_edges.size(); e++) {
            const BSPEdge& edge = complex->edges[f_edges[e]];
            if (GenericPoint3D::point_in_inner_segment(
                    *vertices[tri[i]], *vertices[edge.vertices[0]], *vertices[edge.vertices[1]], xyz
                )) {
                mask |= (1 << i);
                break;
            }
        }
    }
    if (mask == 7) return true;
    for (uint32_t i = 0; i < 3; i++) {
        const uint32_t ti0 = (i + 1) % 3;
        const uint32_t ti1 = (i + 2) % 3;
        if ((mask & (1 << ti0)) != 0 && (mask & (1 << ti1)) != 0) {
            continue;
        }
        for (uint32_t e = 0; e < f_edges.size(); e++) {
            const BSPEdge& edge = complex->edges[f_edges[e]];
            const GenericPoint3D* ev1 = vertices[edge.vertices[0]];
            const GenericPoint3D* ev2 = vertices[edge.vertices[1]];
            const GenericPoint3D* fv1 = vertices[tri[ti0]];
            const GenericPoint3D* fv2 = vertices[tri[ti1]];
            if (GenericPoint3D::inner_segments_cross(*ev1, *ev2, *fv1, *fv2, xyz)) {
                return true;
            }
        }
    }
    return false;
}

inline FaceColor get_face_color(const BSPComplex* complex, const BSPFace& face, const int xyz) {
    const uint32_t* last = complex->edges[face.edges.back()].vertices.data();
    const uint32_t* first = complex->edges[face.edges[0]].vertices.data();
    uint32_t vid = first[0] == last[0] || first[0] == last[1] ? first[0] : first[1];
    const auto& vertices = complex->vertices;
    for (const uint32_t eid : face.edges) {
        const uint32_t* e_verts = complex->edges[eid].vertices.data();
        vid = e_verts[0] == vid ? e_verts[1] : e_verts[0];
        uint32_t out_from_all = 0;
        for (const uint32_t cid : face.coplanar_constraints) {
            const uint32_t* tri = &complex->constraints->triangles[cid * 3];
            if (vid == tri[0] || vid == tri[1] || vid == tri[2]) {
                break;
            }
            const int lpt = localized_point_on_triangle(
                *vertices[vid], *vertices[tri[0]], *vertices[tri[1]], *vertices[tri[2]], xyz
            );
            if (lpt == 2) {
                return FaceColor::BLACK;
            } else if (lpt == 1) {
                break;
            } else {
                out_from_all += 1;
            }
        }
        if (out_from_all == face.coplanar_constraints.size()) {
            return FaceColor::WHITE;
        }
    }

    for (const uint32_t cid : face.coplanar_constraints) {
        const uint32_t* tri = &complex->constraints->triangles[cid * 3];
        if (coplanar_constraint_inner_intersects_face(complex, face.edges, tri, xyz)) {
            return FaceColor::BLACK;
        }
    }

    return FaceColor::WHITE;
}

void BSPComplex::decide_color() {
    for (uint32_t fid = 0; fid < faces.size(); fid++) {
        BSPFace& face = faces[fid];
        if (face.color != FaceColor::GRAY) {
            continue;
        }
        int xyz = face_normal_dominant_component(this, face);
        double p[3];
        face_baricenter_approx(this, face, p);
        ExplicitPoint3D face_center(p[0], p[1], p[2]);
        if (is_baricenter_on_face(this, face, face_center, xyz)) {
            bool found = false;
            for (const uint32_t cid : face.coplanar_constraints) {
                const uint32_t* tri = &constraints->triangles[cid * 3];
                if (GenericPoint3D::point_in_triangle(
                        face_center, *vertices[tri[0]], *vertices[tri[1]], *vertices[tri[2]], xyz
                    )) {
                    found = true;
                    face.color = FaceColor::BLACK;
                    break;
                }
            }
            if (!found) {
                face.color = FaceColor::WHITE;
            }
        } else {
            face.color = get_face_color(this, face, xyz);
        }
    }
}

inline bool is_first_cell_below_face(BSPComplex* complex, const BSPFace& face) {
    const uint32_t* tri = &complex->constraints->triangles[face.coplanar_constraints[0] * 3];
    const GenericPoint3D* pv1 = complex->vertices[tri[0]];
    const GenericPoint3D* pv2 = complex->vertices[tri[1]];
    const GenericPoint3D* pv3 = complex->vertices[tri[2]];
    const BSPCell& cell = complex->cells[face.cells[0]];
    std::vector<uint32_t> cell_verts, _cell_edges;
    find_cell_verts_and_edges(complex, cell, cell_verts, _cell_edges);
    for (const uint32_t eid : face.edges) {
        complex->vert_visit[complex->edges[eid].vertices[0]] = 1;
        complex->vert_visit[complex->edges[eid].vertices[1]] = 1;
    }
    for (const uint32_t vid : cell_verts) {
        if (complex->vert_visit[vid] == 0) {
            const GenericPoint3D* cv = complex->vertices[vid];
            const int o = GenericPoint3D::orient3d(*cv, *pv1, *pv2, *pv3);
            if (o != 0) {
                for (const uint32_t eid : face.edges) {
                    complex->vert_visit[complex->edges[eid].vertices[0]] = 0;
                    complex->vert_visit[complex->edges[eid].vertices[1]] = 0;
                }
                return o < 0;
            }
        }
    }
    return false;
}

void BSPComplex::complex_partition() {
    // for (uint32_t i = 0; i < vertices.size(); i++) {
    //     vert_visit[i] = 0;
    // }
    std::vector<double> approx_coords(vertices.size() * 3);
    for (uint32_t i = 0; i < vertices.size(); i++) {
        vertices[i]->to_double_approx(&approx_coords[i * 3]);
    }
    std::vector<uint32_t> face_verts;
    const auto approx_face_area = [&approx_coords, &face_verts, this](const BSPFace& face) {
        using Vec3 = Eigen::Map<const Eigen::Vector3d>;
        face_verts.clear();
        face_vertices(this, face, face_verts);
        const auto tv0 = Vec3(&approx_coords[face_verts[0] * 3]);

        double area = 0.0;
        for (uint32_t i = 2; i < face_verts.size(); i++) {
            const auto tv1 = Vec3(&approx_coords[face_verts[i - 1] * 3]);
            const auto tv2 = Vec3(&approx_coords[face_verts[i] * 3]);
            const auto v1 = tv1 - tv0;
            const auto v2 = tv2 - tv0;
            // the polygon is convex, no need to judge sign
            area += v1.cross(v2).norm();
        }
        return area;
    };
    std::vector<double> face_areas(faces.size(), 0.0);
    for (uint32_t i = 0; i < faces.size(); i++) {
        face_areas[i] = approx_face_area(faces[i]);
    }
    const double total_area = std::accumulate(face_areas.begin(), face_areas.end(), 0.0);
    for (uint32_t i = 0; i < faces.size(); i++) {
        face_areas[i] /= total_area;
    }

    std::vector<uint8_t> evs(edges.size(), 0);
    for (const BSPFace& f : faces) {
        if (f.color == FaceColor::BLACK) {
            for (const uint32_t eid : f.edges) {
                if (evs[eid] < 2) {
                    evs[eid] += 1;
                }
            }
        }
    }

    // vvs == 1 if vertex is on boundary of skin
    std::vector<uint8_t> vvs(vertices.size(), 0);
    for (uint32_t i = 0; i < edges.size(); i++) {
        if (evs[i] == 1) {
            const BSPEdge& e = edges[i];
            vvs[e.vertices[0]] = vvs[e.vertices[1]] = 1;
        }
    }
    std::vector<double> cell_costs_external(cells.size() + 1, 0.0);
    std::vector<double> cell_costs_internal(cells.size() + 1, 0.0);
    for (uint32_t i = 0; i < faces.size(); i++) {
        const BSPFace& f = faces[i];
        uint32_t cell1 = f.cells[0];
        uint32_t cell2 = f.cells[1];
        if (cell2 == TriFace::INVALID) {
            cell2 = static_cast<uint32_t>(cells.size());
        }
        if (f.color == FaceColor::BLACK) {
            if (is_first_cell_below_face(this, f)) {
                cell_costs_external[cell1] += face_areas[i];
                cell_costs_internal[cell2] += face_areas[i];
            } else {
                cell_costs_external[cell2] += face_areas[i];
                cell_costs_internal[cell1] += face_areas[i];
            }
        }
    }

    cell_costs_internal[cells.size()] = 1.0;
    GraphCut g(
        static_cast<uint32_t>(cell_costs_external.size()), cell_costs_external.data(), cell_costs_internal.data()
    );
    constexpr double w = 1.1;
    for (uint32_t i = 0; i < faces.size(); i++) {
        BSPFace& face = faces[i];
        if (face.color == FaceColor::BLACK) {
            continue;
        }
        bool is_boundary = true;
        for (const uint32_t eid : face.edges) {
            if (vvs[edges[eid].vertices[0]] == 0 || vvs[edges[eid].vertices[1]] == 0) {
                is_boundary = false;
                break;
            }
        }
        uint32_t cell1 = face.cells[0];
        uint32_t cell2 = face.cells[1];
        if (cell2 == TriFace::INVALID) {
            cell2 = static_cast<uint32_t>(cells.size());
        }
        if (is_boundary) {
            g.add_edge(cell1, cell2, face_areas[i] * w, face_areas[i] * w);
        } else {
            g.add_edge(cell1, cell2, face_areas[i] * w, face_areas[i] * w);
        }
    }
    g.max_flow();
    for (uint32_t i = 0; i < cells.size(); i++) {
        cells[i].place = !g.is_sink[i];
    }
}

struct EdgeGroup {
    std::vector<int> edges;
    uint32_t index{TriFace::INVALID};
    EdgeGroup() {}
    EdgeGroup(std::vector<int>&& vec, const uint32_t i) : edges{std::move(vec)}, index{i} {}
};

static inline uint32_t edge_index(const int hid) { return static_cast<uint32_t>(std::abs(hid)) - 1; }

static inline std::vector<int>
make_edge_group(const BSPComplex* complex, std::vector<int>&& outlines, std::vector<uint32_t>& edge_in_group) {
    const auto on_same_line = [&edge_in_group, complex](const int ha, const int hb) {
        // edge id
        const auto ea = edge_index(ha);
        const auto eb = edge_index(hb);
        // group id
        const auto ga = edge_in_group[ea];
        const auto gb = edge_in_group[eb];

        const auto& vertices = complex->vertices;
        const auto& edges = complex->edges;
        if (ga != TriFace::INVALID && gb != TriFace::INVALID) {
            return ga == gb;
        } else {
            // parents
            const auto pa = edges[ea].parent;
            const auto pb = edges[eb].parent;
            if (pa == pb) {
                return true;
            } else {
                const auto v1 = edges[ea].vertices[0];
                const auto v2 = edges[ea].vertices[1];
                const auto v3 = hb < 0 ? edges[eb].vertices[0] : edges[eb].vertices[1];
                return !GenericPoint3D::mis_alignment(*vertices[v1], *vertices[v2], *vertices[v3]);
            }
        }
    };

    std::vector<std::vector<int>> face_edge_groups;
    std::vector<int> group;
    for (const auto hid : outlines) {
        if (group.empty() || on_same_line(group.back(), hid)) {
            group.emplace_back(hid);
        } else {
            face_edge_groups.emplace_back(std::move(group));
            group = {hid};
        }
    }
    if (on_same_line(group.back(), face_edge_groups[0].front())) {
        std::copy(group.begin(), group.end(), std::back_inserter(face_edge_groups[0]));
    } else {
        face_edge_groups.emplace_back(std::move(group));
    }
}

static inline void merge_face_edges(
    BSPComplex* complex, const std::vector<uint32_t>& face_groups,
    const std::vector<std::pair<uint32_t, bool>>& kept_faces,
    const std::vector<std::vector<std::pair<uint32_t, bool>>>& ef_map,
    std::vector<std::vector<std::pair<uint32_t, bool>>> fe_map, std::vector<EdgeGroup>& edge_groups,
    std::vector<uint32_t>& edge_in_group
) {
    std::unordered_map<uint32_t, int> edge_map;
    for (const uint32_t fid_ind : face_groups) {
        std::uint32_t fid;
        bool reversed;
        std::tie(fid, reversed) = kept_faces[fid_ind];
        for (const auto& pair : fe_map[fid]) {
            const std::uint32_t eid = pair.first;
            const bool e_rev = pair.second ^ reversed;
            const auto it = edge_map.find(eid);
            if (it != edge_map.end()) {
                it->second += e_rev ? -1 : 1;
            } else {
                edge_map.emplace(eid, e_rev ? -1 : 1);
            }
        }
    }

    std::vector<int> half_edges;
    std::unordered_map<uint32_t, int> vert_out_edge_map;
    const auto& edges = complex->edges;
    for (const auto& pair : edge_map) {
        if (pair.second == 0) {
            continue;
        }
        const auto eid = pair.first;
        const auto hid = static_cast<int>(eid + 1) * (pair.second < 0 ? -1 : 1);
        half_edges.emplace_back(hid);
        const auto& edge = edges[eid];
        const auto va = pair.second ? edge.vertices[1] : edge.vertices[0];
        vert_out_edge_map.emplace(va, hid);
    }

    auto& edge_visit = complex->edge_visit;
    std::vector<std::vector<int>> face_edge_groups;
    const auto end_vertex = [&edges](const int hid) {
        const auto eid = edge_index(hid);
        const auto& edge = edges[eid];
        return hid < 0 ? edge.vertices[0] : edge.vertices[1];
    };
    for (const auto hid : half_edges) {
        const auto eid = edge_index(hid);
        if (edge_visit[eid] != 0) {
            continue;
        }
        edge_visit[eid] = 1;
        std::vector<int> outlines{hid};
        do {
            const auto vb = end_vertex(outlines.back());
            const auto next_hid = vert_out_edge_map[vb];
            const auto next_eid = edge_index(next_hid);
            if (edge_visit[eid] != 0) {
                break;
            } else {
                outlines.emplace_back(next_hid);
                edge_visit[next_eid] = 1;
            }
        } while (true);
    }
}

static inline void merge_faces(
    BSPComplex* complex, const uint32_t* constraint_parents, const std::vector<std::pair<uint32_t, bool>>& kept_faces,
    const std::vector<std::vector<std::pair<uint32_t, bool>>>& ef_map
) {
    std::unordered_map<uint32_t, uint32_t> face_map;
    face_map.reserve(kept_faces.size());
    for (const auto& pair : kept_faces) {
        face_map.emplace(pair.first, static_cast<uint32_t>(face_map.size()));
    }

    auto same_plane = [&constraint_parents](const BSPFace& fa, const BSPFace& fb) -> bool {
        if (fa.parent == fb.parent) {
            return true;
        }
        for (const uint32_t fa_cid : fa.coplanar_constraints) {
            const auto& fb_constrs = fb.coplanar_constraints;
            auto it = std::find_if(
                fb_constrs.begin(), fb_constrs.end(),
                [&constraint_parents, fa_cid](const uint32_t fb_cid) {
                    if (constraint_parents[fb_cid] == constraint_parents[fa_cid]) {
                        return true;
                    } else {
                        return false;
                    }
                }
            );
            if (it != fb_constrs.end()) {
                return true;
            }
        }
        return false;
    };

    DisjointSet ds(static_cast<uint32_t>(kept_faces.size()));
    for (const auto& pa : kept_faces) {
        const uint32_t fa_id = pa.first;
        const auto& fa = complex->faces[fa_id];
        const uint32_t fa_id_ind = face_map[fa_id];
        for (const uint32_t eid : fa.edges) {
            for (const auto& pb : ef_map[eid]) {
                const uint32_t fb_id = pb.first;
                if (fa_id == fb_id) {
                    continue;
                }
                const uint32_t fb_id_ind = face_map[fb_id];
                if (fb_id_ind < fa_id_ind) {
                    continue;
                }
                const auto& fb = complex->faces[fb_id];
                if (same_plane(fa, fb)) {
                    ds.merge(fa_id_ind, fb_id_ind);
                }
            }
        }
    }
    std::unordered_map<uint32_t, std::vector<uint32_t>> face_group_map;
    face_group_map.reserve(ds.n_groups);
    for (uint32_t i = 0; i < kept_faces.size(); i++) {
        const auto gid = ds.find_set(i);
        auto it = face_group_map.find(gid);
        if (it != face_group_map.end()) {
            it->second.emplace_back(i);
        } else {
            face_group_map.emplace(gid, std::vector<uint32_t>{i});
        }
    }
}

void BSPComplex::extract_skin(
    const uint32_t* constraint_parents, std::vector<double>& out_points, std::vector<uint32_t>& out_faces,
    std::vector<double>& axes, std::vector<uint32_t>& seperator
) {
    std::vector<uint32_t> kept;                                               // kept faces
    std::vector<std::vector<std::pair<uint32_t, bool>>> ef_map(edges.size()); // pair: <fid, reversed>
    std::vector<std::vector<std::pair<uint32_t, bool>>> fe_map(faces.size()); // pair: <eid, reversed>
    for (uint32_t i = 0; i < faces.size(); i++) {
        const BSPFace& face = faces[i];
        if (face.cells[1] == TriFace::INVALID ? cells[face.cells[0]].place == 0
                                              : cells[face.cells[0]].place == cells[face.cells[1]].place) {
            continue;
        }
        kept.emplace_back(i);
        const uint32_t* last = edges[face.edges.back()].vertices.data();
        const uint32_t* first = edges[face.edges[0]].vertices.data();
        uint32_t pvid = (first[0] == last[0] || first[0] == last[1]) ? first[0] : first[1];
        for (uint32_t j = 0; j < face.edges.size(); j++) {
            const uint32_t eid = face.edges[j];
            if (edge_visit[eid] == 0) {
                edge_visit[eid] = 1;
                // vert_visit[edges[eid].vertices[0]] = 1;
                // vert_visit[edges[eid].vertices[1]] = 1;
            }
            if (edges[eid].vertices[0] == pvid) {
                ef_map[eid].emplace_back(i, false);
                fe_map[i].emplace_back(eid, false);
                pvid = edges[eid].vertices[1];
            } else {
                ef_map[eid].emplace_back(i, true);
                fe_map[i].emplace_back(eid, true);
                pvid = edges[eid].vertices[0];
            }
        }
    }

    // const uint32_t n_points = static_cast<uint32_t>(std::count(vert_visit.begin(), vert_visit.end(), 1));
    // out_points.resize(n_points * 3);
    // std::vector<uint32_t> pmap(vertices.size());
    // uint32_t count = 0;
    // for (uint32_t i = 0; i < vertices.size(); i++) {
    //     if (vert_visit[i] == 0) {
    //         continue;
    //     }
    //     double* p = &out_points[count * 3];
    //     pmap[i] = count++;
    //     vertices[i]->to_double(p);
    // }
    // std::fill(vert_visit.begin(), vert_visit.end(), 0);
    // std::fill(edge_visit.begin(), edge_visit.end(), 0);

    // out_faces.reserve(kept.size() * 3);
    // axes.reserve(kept.size() * 6);
    // seperator.reserve(kept.size() + 1);
    // seperator.emplace_back(0);
    std::vector<bool> face_visit(faces.size(), false);
    std::vector<uint32_t> face_verts;
    std::vector<std::pair<uint32_t, bool>> oriented_faces;
    oriented_faces.reserve(kept.size());
    for (const uint32_t fid : kept) {
        if (face_visit[fid]) {
            continue;
        }
        face_verts.clear();
        const auto& face = faces[fid];
        face_vertices(this, face, face_verts);
        const int xyz = GenericPoint3D::max_component_at_triangle_normal(
            vertices[face.mesh_vertices[0]]->to_explicit().ptr(), vertices[face.mesh_vertices[1]]->to_explicit().ptr(),
            vertices[face.mesh_vertices[2]]->to_explicit().ptr()
        );
        const GenericPoint3D* p0 = vertices[face_verts[0]];
        const GenericPoint3D* p1 = vertices[face_verts[1]];
        const GenericPoint3D* p2 = nullptr;
        for (uint32_t i = 2; i < face_verts.size(); i++) {
            const GenericPoint3D* p = vertices[face_verts[i]];
            if (GenericPoint3D::orient2d(*p0, *p1, *p, xyz) != 0) {
                p2 = p;
                break;
            }
        }

        for (const uint32_t vid : face_verts) {
            vert_visit[vid] = 1;
        }
        const uint32_t inner_cid = cells[face.cells[0]].place == 1 ? face.cells[0] : face.cells[1];
        int ori = 0;
        for (const uint32_t f : cells[inner_cid].faces) {
            if (f == fid) {
                continue;
            }
            std::vector<uint32_t> tmp_face_verts;
            face_vertices(this, faces[f], tmp_face_verts);
            for (const uint32_t vid : tmp_face_verts) {
                if (vert_visit[vid] == 1) {
                    continue;
                }
                ori = GenericPoint3D::orient3d(*p0, *p1, *p2, *vertices[vid]);
                if (ori != 0) {
                    break;
                }
            }
            if (ori != 0) {
                break;
            }
        }
        for (const uint32_t vid : face_verts) {
            vert_visit[vid] = 0;
        }
        std::queue<std::pair<uint32_t, bool>> queue; // <face, reversed>

        const auto push_queue = [&queue, &face_visit, &oriented_faces](const uint32_t fi, bool reversed) {
            queue.emplace(fi, reversed);
            face_visit[fi] = true;
            oriented_faces.emplace_back(fi, reversed);
            /* seperator.emplace_back(seperator.back() + static_cast<uint32_t>(face_verts.size()));
            for (const uint32_t vid : face_verts) {
                out_faces.emplace_back(pmap[vid]);
            }
            axes.resize(axes.size() + 6);
            const BSPFace& f = this->faces[fi];
            using Vec3 = Eigen::Map<Eigen::Vector3d>;
            using CVec3 = const Eigen::Map<const Eigen::Vector3d>;
            CVec3 v0(this->vertices[f.mesh_vertices[0]]->to_explicit().ptr());
            CVec3 v1(this->vertices[f.mesh_vertices[1]]->to_explicit().ptr());
            CVec3 v2(this->vertices[f.mesh_vertices[2]]->to_explicit().ptr());
            Vec3 x_axis(&axes[axes.size() - 6]);
            x_axis = (v1 - v0).normalized();
            const auto z_axis = x_axis.cross((v2 - v0)).normalized();
            Vec3 y_axis(&axes[axes.size() - 3]);
            y_axis = z_axis.cross(x_axis).eval();*/
        };
        // assert( ori != 0);
        push_queue(fid, ori < 0);
        while (!queue.empty()) {
            uint32_t fi;
            bool f_rev;
            std::tie(fi, f_rev) = queue.front();
            queue.pop();
            for (const auto& e_pair : fe_map[fi]) {
                uint32_t eid = e_pair.first;
                if (ef_map[eid].size() != 2) {
                    continue;
                }
                const bool e_rev = f_rev ^ e_pair.second;
                uint32_t nfi;
                bool nrev;
                for (const auto& f_pair : ef_map[eid]) {
                    std::tie(nfi, nrev) = f_pair;
                    if (face_visit[nfi]) {
                        continue;
                    }
                    push_queue(nfi, e_rev == nrev);
                }
            }
        }
    }
    std::fill(edge_visit.begin(), edge_visit.end(), 0);
    merge_faces(this, constraint_parents, oriented_faces, ef_map);
}
