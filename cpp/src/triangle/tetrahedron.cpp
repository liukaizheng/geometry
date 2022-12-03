#include <triangle/tetrahedron.h>

#include <Eigen/Dense>

#include <algorithm>
#include <numeric>
#include <random>
#include <utility>

// "Compact Hilbert Indices", Technical Report CS-2006-07
constexpr std::array<std::array<std::array<uint32_t, 8>, 3>, 8> gray_code() noexcept {
    std::array<std::array<std::array<uint32_t, 8>, 3>, 8> transgc;
    uint32_t gc[8];
    constexpr uint32_t n{3}, N{8}, mask{7};

    for (uint8_t i = 0; i < N; i++) {
        gc[i] = i ^ (i >> 1);
    }

    for (uint32_t e = 0; e < N; e++) {
        for (uint32_t d = 0; d < n; d++) {
            const uint32_t f = e ^ (1 << d);
            const uint32_t travel_bit = e ^ f;
            for (uint32_t i = 0; i < N; i++) {
                const uint32_t k = gc[i] * (travel_bit << 1);
                const uint32_t g = (k | (k / N)) & mask;
                transgc[e][d][i] = g ^ e;
            }
        }
    }
    return transgc;
}

constexpr std::array<std::array<std::array<uint32_t, 8>, 3>, 8> transgc{gray_code()};

constexpr std::array<uint32_t, 8> trailing_set_bits_mod3() noexcept {
    std::array<uint32_t, 8> tsb1mod3;
    tsb1mod3[0] = 0;
    constexpr int n{3}, N{8};
    for (uint8_t i = 1; i < N; i++) {
        int v = ~i;             // Count the 0s.
        v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
        uint32_t c;
        for (c = 0; v; c++) {
            v >>= 1;
        }
        tsb1mod3[i] = c % n;
    }
    return tsb1mod3;
}
constexpr std::array<uint32_t, 8> tsb1mod3{trailing_set_bits_mod3()};

struct SortOption {
    uint32_t threshold = 64;
    uint32_t hilbert_order = 52;
    uint32_t hilbert_limit = 8;
    double ratio = 0.125;
};

inline const double* point(const double* points, const uint32_t index) { return &points[index * 3]; }

inline uint32_t hilbert_split(
    const double* points, uint32_t* indices, const uint32_t n_points, const uint32_t gc0, const uint32_t gc1,
    const double* bbox
) {
    if (n_points == 1) {
        return 0;
    }
    const uint32_t axis = (gc0 ^ gc1) >> 1;
    const double split = (bbox[axis] + bbox[axis + 3]) * .5;

    const int d = ((gc0 & (1 << axis)) == 0) ? 1 : -1;
    uint32_t i = 0, j = 0;
    if (d > 0) {
        while (true) {
            for (; i < n_points; i++) {
                if (point(points, indices[i])[axis] >= split) break;
            }
            for (; j < n_points; j++) {
                if (point(points, indices[n_points - 1 - j])[axis] < split) break;
            }
            if (i + j == n_points) break;
            std::swap(indices[i], indices[n_points - 1 - j]);
        }
    } else {
        while (true) {
            for (; i < n_points; i++) {
                if (point(points, indices[i])[axis] <= split) break;
            }
            for (; j < n_points; j++) {
                if (point(points, indices[n_points - 1 - j])[axis] > split) break;
            }
            if (i + j == n_points) break;
            std::swap(indices[i], indices[n_points - 1 - j]);
        }
    }

    return i;
}

static void hilbert_sort(
    const double* points, uint32_t* indices, const uint32_t n_points, const double* bbox, const SortOption* option,
    const uint32_t e, const uint32_t d, const uint32_t depth
) {
    uint32_t p[9];
    p[0] = 0;
    p[8] = n_points;
    p[4] = hilbert_split(points, indices, p[8], transgc[e][d][3], transgc[e][d][4], bbox);
    p[2] = hilbert_split(points, indices, p[4], transgc[e][d][1], transgc[e][d][2], bbox);
    p[1] = hilbert_split(points, indices, p[2], transgc[e][d][0], transgc[e][d][1], bbox);
    p[3] = hilbert_split(points, &indices[p[2]], p[4] - p[2], transgc[e][d][2], transgc[e][d][3], bbox) + p[2];
    p[6] = hilbert_split(points, &indices[p[4]], p[8] - p[4], transgc[e][d][5], transgc[e][d][6], bbox) + p[4];
    p[5] = hilbert_split(points, &indices[p[4]], p[6] - p[4], transgc[e][d][4], transgc[e][d][5], bbox) + p[4];
    p[7] = hilbert_split(points, &indices[p[6]], p[8] - p[6], transgc[e][d][6], transgc[e][d][7], bbox) + p[6];
    if (option->hilbert_order > 0) {
        if (depth + 1 == option->hilbert_order) {
            return;
        }
    }

    constexpr uint32_t mask {7}, n {3};
    for (uint32_t w = 0; w < 8; w++) {
        if ((p[w + 1] - p[w]) > option->hilbert_limit) {
            uint32_t e_w, k;
            if (w == 0) {
                e_w = 0;
            } else {
                //   calculate e(w) = gc(2 * floor((w - 1) / 2)).
                k = 2 * ((w - 1) / 2);
                e_w = k ^ (k >> 1); // = gc(k).
            }
            k = e_w;
            e_w = ((k << (d + 1)) & mask) | ((k >> (n - d - 1)) & mask);
            const uint32_t ei = e ^ e_w;
            uint32_t d_w;
            if (w == 0) {
                d_w = 0;
            } else {
                d_w = ((w % 2) == 0) ? tsb1mod3[w - 1] : tsb1mod3[w];
            }
            const uint32_t di = (d + d_w + 1) % n;
            // Calculate the bounding box of the sub-box.
            double sbox[6];
            if (transgc[e][d][w] & 1) { // x-axis
                sbox[0] = (bbox[0] + bbox[3]) * 0.5;
                sbox[3] = bbox[3];
            } else {
                sbox[0] = bbox[0];
                sbox[3] = (bbox[0] + bbox[3]) * 0.5;
            }
            if (transgc[e][d][w] & 2) { // y-axis
                sbox[1] = (bbox[1] + bbox[4]) * 0.5;
                sbox[4] = bbox[4];
            } else {
                sbox[1] = bbox[1];
                sbox[4] = (bbox[1] + bbox[4]) * 0.5;
            }
            if (transgc[e][d][w] & 4) { // z-axis
                sbox[2] = (bbox[2] + bbox[5]) * 0.5;
                sbox[5] = bbox[5];
            } else {
                sbox[2] = bbox[2];
                sbox[5] = (bbox[2] + bbox[5]) * 0.5;
            }
            hilbert_sort(points, &indices[p[w]], p[w + 1] - p[w], sbox, option, ei, di, depth + 1);
        }
    }
}

static void brio_multiscale_sort(
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


Tetrahedrons Tetrahedrons::tetrahedralize(const double* points, const uint32_t n_points, const double epsilon) {
    Eigen::Map<const Eigen::Matrix<double, -1, 3, Eigen::RowMajor>> matrix(points, n_points, 3);
    const auto min_corner = matrix.colwise().minCoeff().eval();
    const auto max_corner = matrix.colwise().maxCoeff().eval();
    double bbox[6]{min_corner[0], min_corner[1], min_corner[2], max_corner[0], max_corner[1], max_corner[2]};
    std::vector<uint32_t> sorted_pt_inds(n_points);
    std::iota(sorted_pt_inds.begin(), sorted_pt_inds.end(), 0);
    std::shuffle(sorted_pt_inds.begin(), sorted_pt_inds.end(), std::default_random_engine(2));
    SortOption option;
    uint32_t n_group;
    brio_multiscale_sort(points, sorted_pt_inds.data(), n_points, bbox, &option, n_group);

    Tetrahedrons tets;
    return tets;
}
