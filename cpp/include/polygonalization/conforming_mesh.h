#pragma once

#include <triangle/tetrahedron.h>

#include <vector>

struct Constraints {
    std::vector<uint32_t> triangles;
    const uint32_t n_triangles; // the number of origin triangles
    
    Constraints(const uint32_t* indices, const uint32_t n): n_triangles{n} {
        const uint32_t len = n_triangles * 3;
        triangles.resize(len);
        for (uint32_t i = 0; i < len; i++) {
           triangles[i] = indices[i];
        }
    }
};

void place_virtual_constraints(const TetMesh& mesh, Constraints& constraints);
void insert_constraints(const TetMesh& mesh, const Constraints& constraints, std::vector<std::vector<uint32_t>>* tet_map);
