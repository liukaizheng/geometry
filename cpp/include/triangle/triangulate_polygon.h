#pragma once
#include <vector>

#ifdef EMSCRIPTEN
extern "C" {
#endif

uint32_t triangulate_polygon_soup(
    const double* points, const uint32_t* edge_data, const double* axis_data, const uint32_t* seperator,
    const uint32_t n_polygon, uint32_t** triangles, uint32_t** triangle_parents
); 

#ifdef EMSCRIPTEN
}

#endif
