#pragma once

#include <array>
#include <vector>

struct Tet {
    std::array<uint32_t, 4> data;
};

struct Tetrahedrons {
    std::vector<Tet> tets;
    
    static Tetrahedrons tetrahedralize(const double* points, const uint32_t n_points, const double epsilon);
};
