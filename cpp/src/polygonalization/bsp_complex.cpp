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
                for (const uint32_t ev : {edge.vertices[0], edge.vertices[1]}) {
                    if (vert_visit[ev] == 0) {
                        vert_visit[ev] = 1;
                        cell_verts.emplace_back(ev);
                    }
                }
            }
        }
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
        if (consecutive_edges(new_edge.vertices, complex->edges[face_edges[i]].vertices)) {
            const uint32_t pos = (i + 1) % n_edges;
            if (consecutive_edges(new_edge.vertices, complex->edges[face_edges[pos]].vertices)) {
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
        }
    }
}
