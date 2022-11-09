#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include <vector>

extern "C" {
void exactinit();
uint32_t* triangulate_polygon_soup(
    const double* points, const uint32_t* edge_data, const double* axis_data, const uint32_t* seperator,
    const uint32_t n_polygon, uint32_t* n_triangles
);
}

TEST_CASE("test: triangulate a cube") {
    std::vector<double> point_data {
        0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1,
    };
    std::vector<uint32_t> indices {
        3, 2, 2, 1, 1, 0, 0, 3, 4, 5, 5, 6, 6, 7, 7, 4, 0, 1, 1, 5, 5, 4, 4, 0,
        1, 2, 2, 6, 6, 5, 5, 1, 3, 7, 7, 6, 6, 2, 2, 3, 4, 7, 7, 3, 3, 0, 0, 4,
    };
    std::vector<uint32_t> seperator{0, 8, 16, 24, 32, 40, 48};
    std::vector<double> axes{
        0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1,  0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
    };
    uint32_t n_triangles = 0;
    const auto result = triangulate_polygon_soup(
        point_data.data(), indices.data(), axes.data(), seperator.data(), static_cast<uint32_t>(seperator.size() - 1),
        &n_triangles
    );
    delete[] result;
    CHECK(n_triangles ==  12);
}
