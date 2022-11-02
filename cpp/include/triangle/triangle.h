#pragma once

#include <vector>

#ifdef EMSCRIPTEN
extern "C" {
#endif

uint32_t* triangulate(const double* points, const uint32_t n_points, const uint32_t* segments, const uint32_t n_segments, uint32_t* n_triangles); 

#ifdef EMSCRIPTEN
}
#endif

