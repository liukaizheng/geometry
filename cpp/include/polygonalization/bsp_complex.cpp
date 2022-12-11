#include <polygonalization/bsp_complex.h>

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
            face.color = FaceColor::BLACK;
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
}
