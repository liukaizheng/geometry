#include <triangle/tetrahedron.h>

#include <algorithm>
#include <numeric>
#include <random>

// "Compact Hilbert Indices", Technical Report CS-2006-07
constexpr std::array<std::array<std::array<int, 8>, 3>, 8> gray_code() noexcept {
    std::array<std::array<std::array<int, 8>, 3>, 8> transgc;
    int gc[8];
    constexpr int n{3}, N{8}, mask{7};

    for (int i = 0; i < N; i++) {
        gc[i] = i ^ (i >> 1);
    }

    for (uint8_t e = 0; e < N; e++) {
        for (uint8_t d = 0; d < n; d++) {
            const int f = e ^ (1 << d);
            const int travel_bit = e ^ f;
            for (uint8_t i = 0; i < N; i++) {
                const int k = gc[i] * (travel_bit << 1);
                const int g = (k | (k / N)) & mask;
                transgc[e][d][i] = g ^ e;
            }
        }
    }
    return transgc;
}

constexpr std::array<std::array<std::array<int, 8>, 3>, 8> transgc{gray_code()};

struct SortOption {
    uint32_t threshold = 64;
    uint32_t hilbert_order = 52;
    uint32_t hilbert_limit = 8;
    double ratio = 0.125;
};
uint32_t hilbert_split(
    const double* points, uint32_t* indices, const uint32_t n_points, const uint32_t gc0, const uint32_t gc1,
    const double* bbox
) {
    if (n_points == 1) {
        return 0;
    }
    const uint32_t axis = (gc0 ^ gc1) >> 1;
    const double split = (bbox[axis] + bbox[axis + 3]) * .5;

    const int d = ((gc0 & (1 << axis)) == 0) ? 1 : -1;
    if (d > 0) {
        uint32_t i = 0, j = n_points - 1;
        do {
            for (uint32_t i = 0; i < n_points; i++) {
                if (points[indices[i] * 3 + axis] >= split) break;
            }
            for (; j >= 0; j--) {
                if (vertexarray[j][axis] < split) break;
            }
            // Is the partition finished?
            if (i == (j + 1)) break;
            // Swap i-th and j-th vertices.
            swapvert = vertexarray[i];
            vertexarray[i] = vertexarray[j];
            vertexarray[j] = swapvert;
            // Continue patitioning the array;
        } while (true);
    } else {
        do {
            for (; i < arraysize; i++) {
                if (vertexarray[i][axis] <= split) break;
            }
            for (; j >= 0; j--) {
                if (vertexarray[j][axis] > split) break;
            }
            // Is the partition finished?
            if (i == (j + 1)) break;
            // Swap i-th and j-th vertices.
            swapvert = vertexarray[i];
            vertexarray[i] = vertexarray[j];
            vertexarray[j] = swapvert;
            // Continue patitioning the array;
        } while (true);
    }

    return i;
}

void hilbert_sort(
    const double* points, uint32_t* indices, const uint32_t n_points, const double* bbox, const SortOption* option,
    const uint32_t e, const uint32_t d, const uint32_t depth
) {}

void brio_multiscale_sort(
    const double* points, uint32_t* indices, const uint32_t n_points, const double* bbox, const SortOption* option,
    uint32_t& depth
) {
    uint32_t middle = 0;
    if (n_points >= option->threshold) {
        depth += 1;
        middle = static_cast<uint32_t>(static_cast<double>(n_points) * option->ratio);
        brio_multiscale_sort(points, indices, middle, bbox, option, depth);
    }
    hilbert_sort(points, indices + middle, n_points - middle, bbox, option, 0, 0, 0);
}

constexpr std::array<int, 8> trailing_set_bits_mod3() noexcept {
    std::array<int, 8> tsb1mod3;
    tsb1mod3[0] = 0;
    constexpr int n{3}, N{8};
    for (uint8_t i = 1; i < N; i++) {
        int v = ~i;             // Count the 0s.
        v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
        int c;
        for (c = 0; v; c++) {
            v >>= 1;
        }
        tsb1mod3[i] = c % n;
    }
    return tsb1mod3;
}
constexpr std::array<int, 8> tsb1mod3{trailing_set_bits_mod3()};

Tetrahedrons Tetrahedrons::tetrahedralize(const double* points, const uint32_t n_points, const double epsilon) {
    Tetrahedrons tets;
    std::vector<uint32_t> sorted_pt_inds(n_points);
    std::iota(sorted_pt_inds.begin(), sorted_pt_inds.end(), 0);
    std::shuffle(sorted_pt_inds.begin(), sorted_pt_inds.end(), std::default_random_engine(2));
    return tets;
}
