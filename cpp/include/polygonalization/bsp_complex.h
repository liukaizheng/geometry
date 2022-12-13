#pragma once

#include <polygonalization/conforming_mesh.h>
#include <predicates/generic_point.h>
#include <triangle/tetrahedron.h>

enum class FaceColor { BLACK, WHITE, GRAY };

struct BSPEdge {
    uint32_t mesh_vertices[6]{TriFace::INVALID, TriFace::INVALID, TriFace::INVALID,
                              TriFace::INVALID, TriFace::INVALID, TriFace::INVALID};
    std::array<uint32_t, 2> vertices;
    uint32_t face;

    BSPEdge() {}
    BSPEdge(const uint32_t ev1, const uint32_t ev2) {
        mesh_vertices[0] = vertices[0] = ev1;
        mesh_vertices[1] = vertices[1] = ev2;
    }
    BSPEdge(const uint32_t ev1, const uint32_t ev2, const uint32_t* tri1, const uint32_t* tri2) {
        vertices[0] = ev1;
        vertices[1] = ev2;
        mesh_vertices[0] = tri1[0];
        mesh_vertices[1] = tri1[1];
        mesh_vertices[2] = tri1[2];
        mesh_vertices[3] = tri2[0];
        mesh_vertices[4] = tri2[1];
        mesh_vertices[5] = tri2[2];
    }

    BSPEdge split(const uint32_t vid) {
        BSPEdge e;
        if (mesh_vertices[2] == TriFace::INVALID) {
            e.mesh_vertices[0] = mesh_vertices[0];
            e.mesh_vertices[1] = mesh_vertices[1];
        } else {
            e.mesh_vertices[0] = mesh_vertices[0];
            e.mesh_vertices[1] = mesh_vertices[1];
            e.mesh_vertices[2] = mesh_vertices[2];
            e.mesh_vertices[3] = mesh_vertices[3];
            e.mesh_vertices[4] = mesh_vertices[4];
            e.mesh_vertices[5] = mesh_vertices[5];
        }
        e.vertices[0] = vertices[0];
        e.vertices[1] = vid;
        vertices[0] = vid;
        e.face = face;
        return e;
    }
};

struct BSPFace {
    std::vector<uint32_t> edges;
    uint32_t cells[2];
    std::vector<uint32_t> coplanar_constraints;
    uint32_t mesh_vertices[3];
    FaceColor color;

    BSPFace() {}
    BSPFace(const uint32_t* verts, const uint32_t cell_ind, const uint32_t adj_cell_ind)
        : cells{cell_ind, adj_cell_ind}, mesh_vertices{verts[0], verts[1], verts[1]} {}
};

struct BSPCell {
    std::vector<uint32_t> faces;
    std::vector<uint32_t> constraints;
    
    BSPCell() {}
    
    void remove_face(const uint32_t pos) {
        if (pos + 1 != faces.size()) {
            faces[pos] = faces.back();
        }
        faces.pop_back();
    }
};

struct BSPComplex {
    std::vector<GenericPoint3D*> vertices;
    std::vector<BSPEdge> edges;
    std::vector<BSPFace> faces;
    std::vector<BSPCell> cells;
    const Constraints* constraints;
    std::vector<int> verts_oris;
    std::vector<uint32_t> vert_visit;
    std::vector<uint32_t> edge_visit;

    BSPComplex(
        const TetMesh& mesh, const Constraints* constraints,
        std::array<std::vector<std::vector<uint32_t>>, 5>&& tet_maps
    );
    ~BSPComplex() {
        for (GenericPoint3D* p : vertices) {
            delete p;
        }
    }

    void split_cell(const uint32_t cid);
};
