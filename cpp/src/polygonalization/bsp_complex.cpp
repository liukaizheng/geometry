#include <iterator>
#include <polygonalization/bsp_complex.h>

#include <algorithm>
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
        edges.emplace_back(ev1, ev2);
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
    const TetMesh& mesh, const Constraints* c, std::array<std::vector<std::vector<uint32_t>>, 5>&& tet_maps
)
    : constraints(c) {
    vertices.reserve(mesh.n_points);
    for (uint32_t i = 0; i < mesh.n_points; i++) {
        const double* p = mesh.point(i);
        vertices.emplace_back(new ExplicitPoint3D(p[0], p[1], p[2]));
    }
    verts_oris.resize(vertices.size() << 1, 2);
    std::vector<uint32_t> new_order;
    const uint32_t n_cells = remove_ghost_tets(mesh, new_order);
    cells.resize(n_cells);
    edges.reserve(n_cells + mesh.n_points);
    faces.reserve(cells.size() << 1);

    std::vector<std::unordered_map<uint32_t, uint32_t>> edge_map(mesh.n_points - 1);
    std::vector<uint32_t> tet_edges;
    tet_edges.reserve(6);
    for (uint32_t i = 0; i < mesh.tets.size(); i++) {
        const uint32_t cell_ind = new_order[i];
        if (cell_ind == TriFace::INVALID) {
            continue;
        }
        cells[cell_ind].constraints = tet_maps[4][i];
        constexpr std::array<uint32_t, 12> edgetbl{{3, 5, 4, 1, 5, 2, 2, 4, 0}};
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
                const uint32_t face_ind = static_cast<uint32_t>(faces.size());
                faces.emplace_back(verts, cell_ind, adj_cell_ind);
                cells[cell_ind].faces.emplace_back(face_ind);
                if (adj_cell_ind != TriFace::INVALID) {
                    cells[adj_cell_ind].faces.emplace_back(face_ind);
                }
                const uint32_t* edge_indices = &edgetbl[j * 3];
                for (uint32_t k = 0; k < 3; k++) {
                    const uint32_t eid = tet_edges[edge_indices[k]];
                    faces.back().edges.emplace_back(eid);
                    edges[eid].face = face_ind;
                }
                fill_face_color(i, face_ind, tet_maps[j], this);
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
    std::vector<uint32_t>& edge_visit = complex->edge_visit;
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
        vertices[edge.vertices[0]]->to_explicit(), vertices[edge.vertices[1]]->to_explicit(),
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
        comm[i++] = v1[1];
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

inline std::vector<uint32_t> face_vertices(BSPComplex* complex, const BSPFace& face) {
    auto& vert_visit = complex->vert_visit;
    auto& edges = complex->edges;
    std::vector<uint32_t> result(face.edges.size());
    uint32_t idx = 0;
    for (const uint32_t eid : face.edges) {
        const BSPEdge& edge = edges[eid];
        for (const uint32_t vid : edge.vertices) {
            if (vert_visit[vid] == 0) {
                vert_visit[vid] = 1;
                result[idx++] = vid;
            }
        }
    }
    for (const uint32_t vid : result) {
        vert_visit[vid] = 0;
    }
    return result;
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
    complex->faces[new_fid].edges.assign(face.edges.begin(), face.edges.end());
    face.edges.resize(idx);
    for (const uint32_t eid : face.edges) {
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
        endpts[0], endpts[1], face.mesh_vertices, &complex->constraints->triangles[constr_id * 3]
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
    complex->faces.emplace_back(face.mesh_vertices, face.cells[0], face.cells[1]);
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
            zero_verts[pos++] = 0;
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
        const auto face_verts = face_vertices(complex, face);
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
            if (v_edges[0] != TriFace::INVALID) {
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
    complex->faces.emplace_back(&complex->constraints->triangles[constr_id * 3], cid, new_cid);
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
        const uint32_t* triangle = &complex->constraints->triangles[down_cell.constraints[i]];
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

    for (const uint32_t eid : cell_edges) {
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
        const std::vector<uint32_t> face_verts = face_vertices(this, face);
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
    for (const uint32_t vid : cell_verts) {
        verts_oris[vid] = 2;
    }
    auto reset_oris = [this](const std::vector<uint32_t>& constrs) {
        for (const uint32_t c : constrs) {
            const uint32_t* tri = &constraints->triangles[c * 3];
            verts_oris[tri[0]] = 2;
            verts_oris[tri[1]] = 2;
            verts_oris[tri[2]] = 2;
        }
    };
    reset_oris(faces.back().coplanar_constraints);
    reset_oris(cells[cid].constraints);
    reset_oris(cells[new_cid].constraints);
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
    const BSPComplex* complex, const std::vector<uint32_t>& f_edges, const uint32_t* tri, int xyz
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
    /* todo */
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
                const double* p0 = vertices[tri[0]]->to_explicit().ptr();
                const double* p1 = vertices[tri[1]]->to_explicit().ptr();
                const double* p2 = vertices[tri[2]]->to_explicit().ptr();
                if (GenericPoint3D::point_in_triangle(p, p0, p1, p2)) {
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
