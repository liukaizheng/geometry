#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>

void make_polyhedral_mesh_from_triangles(
    const double* points, const uint32_t n_points, const uint32_t* triangles, const uint32_t n_triangles, const std::vector<std::unordered_map<uint32_t, uint32_t>>& ori_edge_parents,
    std::vector<double>& out_points, std::vector<uint32_t>& out_faces, std::vector<double>& axes,
    std::vector<uint32_t>& seperators
);
