#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

void make_polyhedral_mesh_from_triangles(
    const double* points, const uint32_t n_points, const uint32_t* triangles, const uint32_t n_triangles,
    const std::vector<std::unordered_map<uint32_t, uint32_t>>& ori_edge_parents, const uint32_t* face_parents,
    std::vector<double>& out_points, std::vector<uint32_t>& out_loops, std::vector<uint32_t>& out_loop_separators,
    std::vector<uint32_t>& out_faces, std::vector<uint32_t>& out_face_separators, std::vector<double>& axes
);
