#include <predicates/predicates.h>
#include <triangle/tetrahedron.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
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

    constexpr uint32_t mask{7}, n{3};
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

constexpr uint32_t INVALID{std::numeric_limits<uint32_t>::max()};

constexpr std::array<std::array<uint32_t, 12>, 12> bond_table() noexcept {
    std::array<std::array<uint32_t, 12>, 12> bondtbl;
    for (uint32_t i = 0; i < 12; i++) {
        for (uint32_t j = 0; j < 12; j++) {
            bondtbl[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
        }
    }
    return bondtbl;
}
constexpr std::array<std::array<uint32_t, 12>, 12> bondtbl{bond_table()};

constexpr std::array<uint32_t, 12> esymtbl{{9, 6, 11, 4, 3, 7, 1, 5, 10, 0, 8, 2}};

constexpr std::array<uint32_t, 12> enext_table() noexcept {
    std::array<uint32_t, 12> enexttbl;
    for (uint32_t i = 0; i < 12; i++) {
        enexttbl[i] = (i + 4) % 12;
    }
    return enexttbl;
}
constexpr std::array<uint32_t, 12> enexttbl{enext_table()};

constexpr std::array<uint32_t, 12> eprev_table() noexcept {
    std::array<uint32_t, 12> eprevtbl;
    for (uint32_t i = 0; i < 12; i++) {
        eprevtbl[i] = (i + 8) % 12;
    }
    return eprevtbl;
}
constexpr std::array<uint32_t, 12> eprevtbl{eprev_table()};

constexpr std::array<uint32_t, 12> enext_esym_table() noexcept {
    std::array<uint32_t, 12> enextesymtbl;
    for (uint32_t i = 0; i < 12; i++) {
        enextesymtbl[i] = esymtbl[enexttbl[i]];
    }
    return enextesymtbl;
}
constexpr std::array<uint32_t, 12> enextesymtbl{enext_esym_table()};

constexpr std::array<uint32_t, 12> eprev_esym_table() noexcept {
    std::array<uint32_t, 12> eprevesymtbl;
    for (uint32_t i = 0; i < 12; i++) {
        eprevesymtbl[i] = esymtbl[eprevtbl[i]];
    }
    return eprevesymtbl;
}
constexpr std::array<uint32_t, 12> eprevesymtbl{eprev_esym_table()};

constexpr std::array<uint32_t, 12> orgpivot{{3, 3, 1, 1, 2, 0, 0, 2, 1, 2, 3, 0}};
constexpr std::array<uint32_t, 12> destpivot{{2, 0, 0, 2, 1, 2, 3, 0, 3, 3, 1, 1}};
constexpr std::array<uint32_t, 12> apexpivot{{1, 2, 3, 0, 3, 3, 1, 1, 2, 0, 0, 2}};
constexpr std::array<uint32_t, 12> oppopivot{{0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3}};

constexpr std::array<std::array<uint32_t, 12>, 12> fsym_table() noexcept {
    std::array<std::array<uint32_t, 12>, 12> fsymtbl;
    for (uint32_t i = 0; i < 12; i++) {
        for (uint32_t j = 0; j < 12; j++) {
            fsymtbl[i][j] = (j + 12 - (i & 12)) % 12;
        }
    }
    return fsymtbl;
}
constexpr std::array<std::array<uint32_t, 12>, 12> fsymtbl{fsym_table()};

constexpr std::array<uint32_t, 12> face_pivot1() noexcept {
    std::array<uint32_t, 12> facepivot1;
    for (uint32_t i = 0; i < 12; i++) {
        facepivot1[i] = (esymtbl[i] & 3);
    }
    return facepivot1;
}
constexpr std::array<uint32_t, 12> facepivot1{face_pivot1()};

constexpr std::array<std::array<uint32_t, 12>, 12> face_pivot2() noexcept {
    std::array<std::array<uint32_t, 12>, 12> facepivot2;
    for (uint32_t i = 0; i < 12; i++) {
        for (uint32_t j = 0; j < 12; j++) {
            facepivot2[i][j] = fsymtbl[esymtbl[i]][j];
        }
    }
    return facepivot2;
}
constexpr std::array<std::array<uint32_t, 12>, 12> facepivot2{face_pivot2()};

inline void bond(Tetrahedrons& tets, const TriFace& t1, const TriFace& t2) {
    TriFace& f1 = tets.tets[t1.tet].nei[t1.ver & 3];
    f1.tet = t2.tet;
    f1.ver = bondtbl[t1.ver][t2.ver];

    TriFace& f2 = tets.tets[t2.tet].nei[t2.ver & 3];
    f2.tet = t1.tet;
    f2.ver = bondtbl[t2.ver][t1.ver];
}

// next edge
inline void enext(const TriFace& t1, TriFace& t2) {
    t2.tet = t1.tet;
    t2.ver = enexttbl[t1.ver];
}

inline void enextself(TriFace& t) { t.ver = enexttbl[t.ver]; }

inline void eprev(const TriFace& t1, TriFace& t2) {
    t2.tet = t1.tet;
    t2.ver = eprevtbl[t1.ver];
}

inline void eprevself(TriFace& t) { t.ver = eprevtbl[t.ver]; }

// symmetric edge
inline void esym(const TriFace& t1, TriFace& t2) {
    t2.tet = t1.tet;
    t2.ver = esymtbl[t1.ver];
}

inline void esymself(TriFace& t) { t.ver = esymtbl[t.ver]; }

inline void enextesym(const TriFace& t1, TriFace& t2) {
    t2.tet = t1.tet;
    t2.ver = enextesymtbl[t1.ver];
}

inline void enextesymself(TriFace& t) { t.ver = enextesymtbl[t.ver]; }

inline void eprevesym(const TriFace& t1, TriFace& t2) {
    t2.tet = t1.tet;
    t2.ver = eprevesymtbl[t1.ver];
}

inline void eprevesymself(TriFace& t) { t.ver = eprevesymtbl[t.ver]; }

inline void fsym(const Tetrahedrons& tets, const TriFace& t1, TriFace& t2) {
    const TriFace& nf = tets.tets[t1.tet].nei[t1.ver & 3];
    t2.tet = nf.tet;
    t2.ver = fsymtbl[t1.ver][nf.ver];
}

inline void fsymself(const Tetrahedrons& tets, TriFace& t) {
    const TriFace& nf = tets.tets[t.tet].nei[t.ver & 3];
    t.tet = nf.tet;
    t.ver = fsymtbl[t.ver][nf.ver];
}

inline void fnextself(const Tetrahedrons& tets, TriFace& t) {
    const TriFace& nf = tets.tets[t.tet].nei[facepivot1[t.ver]];
    t.tet = nf.tet;
    t.ver = facepivot2[t.ver][nf.ver];
}

inline void infect(Tetrahedrons& tets, const uint32_t t) { tets.tets[t].mask |= 1; }

inline void uninfect(Tetrahedrons& tets, const uint32_t t) { tets.tets[t].mask &= ~1; }

inline bool infected(const Tetrahedrons& tets, const uint32_t t) { return (tets.tets[t].mask & 1) != 0; }

inline uint32_t org(const Tetrahedrons& tets, const TriFace& f) { return tets.tets[f.tet].data[orgpivot[f.ver]]; }

inline uint32_t dest(const Tetrahedrons& tets, const TriFace& f) { return tets.tets[f.tet].data[destpivot[f.ver]]; }

inline uint32_t apex(const Tetrahedrons& tets, const TriFace& f) { return tets.tets[f.tet].data[destpivot[f.ver]]; }

inline uint32_t oppo(const Tetrahedrons& tets, const TriFace& f) { return tets.tets[f.tet].data[oppopivot[f.ver]]; }

inline double
orient3d(const double* points, const uint32_t pa, const uint32_t pb, const uint32_t pc, const uint32_t pd) {
    return orient3d(&points[pa * 3], &points[pb * 3], &points[pc * 3], &points[pd * 3]);
}

inline void make_tet(
    Tetrahedrons& tets, TriFace& face, const uint32_t pa, const uint32_t pb, const uint32_t pc, const uint32_t pd
) {
    face.tet = static_cast<uint32_t>(tets.tets.size());
    face.ver = 11;
    tets.tets.emplace_back(std::array<uint32_t, 4>{{pa, pb, pc, pd}}, std::array<TriFace, 4>{{{}, {}, {}, {}}});
}

inline TriFace
initial_delaunay(Tetrahedrons& tets, const uint32_t pa, const uint32_t pb, const uint32_t pc, const uint32_t pd) {
    TriFace firsttet, tetopa, tetopb, tetopc, tetopd;

    make_tet(tets, firsttet, pa, pb, pc, pd);
    make_tet(tets, tetopa, pb, pc, pd, tets.n_points);
    make_tet(tets, tetopb, pc, pa, pd, tets.n_points);
    make_tet(tets, tetopc, pa, pb, pd, tets.n_points);
    make_tet(tets, tetopd, pb, pa, pc, tets.n_points);

    // Connect hull tetrahedra to firsttet (at four faces of firsttet).
    TriFace worktet;
    bond(tets, firsttet, tetopd);
    esym(firsttet, worktet);
    bond(tets, worktet, tetopc); // ab
    enextesym(firsttet, worktet);
    bond(tets, worktet, tetopa); // bc
    eprevesym(firsttet, worktet);
    bond(tets, worktet, tetopb); // ca

    // Connect hull tetrahedra together (at six edges of firsttet).
    TriFace worktet1;
    esym(tetopc, worktet);
    esym(tetopd, worktet1);
    bond(tets, worktet, worktet1); // ab
    esym(tetopa, worktet);
    eprevesym(tetopd, worktet1);
    bond(tets, worktet, worktet1); // bc
    esym(tetopb, worktet);
    enextesym(tetopd, worktet1);
    bond(tets, worktet, worktet1); // ca
    eprevesym(tetopc, worktet);
    enextesym(tetopb, worktet1);
    bond(tets, worktet, worktet1); // da
    eprevesym(tetopa, worktet);
    enextesym(tetopc, worktet1);
    bond(tets, worktet, worktet1); // db
    eprevesym(tetopb, worktet);
    enextesym(tetopa, worktet1);
    bond(tets, worktet, worktet1); // dc

    tets.p2t[pa] = 0;
    tets.p2t[pb] = 0;
    tets.p2t[pc] = 0;
    tets.p2t[pd] = 0;

    tets.p2t[tets.n_points] = tetopa.tet;

    return firsttet;
}

enum class LocateResult { OUTSIDE, ONVERTEX, ONEDGE, ONFACE, INTETRAHEDRON };

inline LocateResult locate_dt(const Tetrahedrons& tets, const uint32_t pid, TriFace& searchtet) {
    if (tets.is_hull_tet(searchtet.tet)) {
        searchtet.tet = tets.tets[searchtet.tet].nei[3].tet;
    }

    for (searchtet.ver = 0; searchtet.ver < 4; searchtet.ver++) {
        const double ori =
            orient3d(tets.points, org(tets, searchtet), dest(tets, searchtet), apex(tets, searchtet), pid);
        if (ori < 0.0) break;
    }

    LocateResult loc = LocateResult::OUTSIDE;
    while (true) {
        uint32_t toppo = oppo(tets, searchtet);

        // Check if the vertex is we seek.
        if (toppo == pid) {
            // Adjust the origin of searchtet to be searchpt.
            esymself(searchtet);
            eprevself(searchtet);
            loc = LocateResult::ONVERTEX; // return ONVERTEX;
            break;
        }

        const double oriorg = orient3d(tets.points, dest(tets, searchtet), apex(tets, searchtet), toppo, pid);
        if (oriorg < 0.0) {
            enextesymself(searchtet);
        } else {
            const double oridest = orient3d(tets.points, apex(tets, searchtet), org(tets, searchtet), toppo, pid);
            if (oridest < 0) {
                eprevesymself(searchtet);
            } else {
                const double oriapex = orient3d(tets.points, org(tets, searchtet), dest(tets, searchtet), toppo, pid);
                if (oriapex < 0) {
                    esymself(searchtet);
                } else {
                    // oriorg >= 0, oridest >= 0, oriapex >= 0 ==> found the point.
                    // The point we seek must be on the boundary of or inside this
                    //   tetrahedron. Check for boundary cases first.
                    if (oriorg == 0.0) {
                        // Go to the face opposite to origin.
                        enextesymself(searchtet);
                        if (oridest == 0.0) {
                            eprevself(searchtet); // edge oppo->apex
                            if (oriapex == 0.0) {
                                // oppo is duplicated with p.
                                loc = LocateResult::ONVERTEX; // return ONVERTEX;
                                break;
                            }
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        if (oriapex == 0.0) {
                            enextself(searchtet);       // edge dest->oppo
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        loc = LocateResult::ONFACE; // return ONFACE;
                        break;
                    }
                    if (oridest == 0.0) {
                        // Go to the face opposite to destination.
                        eprevesymself(searchtet);
                        if (oriapex == 0.0) {
                            eprevself(searchtet);       // edge oppo->org
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        loc = LocateResult::ONFACE; // return ONFACE;
                        break;
                    }
                    if (oriapex == 0.0) {
                        // Go to the face opposite to apex
                        esymself(searchtet);
                        loc = LocateResult::ONFACE; // return ONFACE;
                        break;
                    }
                    loc = LocateResult::INTETRAHEDRON;
                    break;
                }
            }
        }

        const TriFace& nf = tets.tets[searchtet.tet].nei[searchtet.ver & 3];
        searchtet.tet = nf.tet;
        searchtet.ver = nf.ver;

        if (tets.is_hull_tet(searchtet.tet)) {
            loc = LocateResult::OUTSIDE;
            break;
        }
    }

    return loc;
}

// Insert a vertex using the Bowyer-Watson algorithm
inline bool insert_vertex_bw(Tetrahedrons& tets, const uint32_t pid, TriFace& searchtet) {
    const LocateResult loc = locate_dt(tets, pid, searchtet);
    std::vector<uint32_t> cave_oldtet_list;
    if (loc == LocateResult::OUTSIDE) {
        infect(tets, searchtet.tet);
        cave_oldtet_list.emplace_back(searchtet.tet);
    } else if (loc == LocateResult::INTETRAHEDRON) {
        infect(tets, searchtet.tet);
        cave_oldtet_list.emplace_back(searchtet.tet);
    } else if (loc == LocateResult::ONFACE) {
        infect(tets, searchtet.tet);
        cave_oldtet_list.emplace_back(searchtet.tet);
        const uint32_t nei_tet = tets.tets[searchtet.tet].nei[searchtet.ver & 3].tet;
        infect(tets, nei_tet);
        cave_oldtet_list.emplace_back(nei_tet);
    } else if (loc == LocateResult::ONEDGE) {
        TriFace spintet = searchtet;
        while (true) {
            infect(tets, spintet.tet);
            cave_oldtet_list.emplace_back(spintet.tet);
            fnextself(tets, spintet);
            if (spintet.tet == searchtet.tet) break;
        }
    } else if (loc == LocateResult::ONVERTEX) {
        // The point already exist. Do nothing and return.
        return false;
    }
    return false;
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

    const double bbox_size = (max_corner - min_corner).norm();
    const double bbox_size2 = bbox_size * bbox_size;
    const double bbox_size3 = bbox_size2 * bbox_size;
    uint32_t i = 1;
    using Vector3d = Eigen::Map<const Eigen::Vector3d>;
    {
        Vector3d p0{point(points, sorted_pt_inds[0])};
        while (i < n_points) {
            Vector3d pi{point(points, sorted_pt_inds[i])};
            if ((p0 - pi).norm() / bbox_size < epsilon) {
                i += 1;
            } else {
                break;
            }
        }
        if (i > 1) {
            if (i == n_points) return {};
            std::swap(sorted_pt_inds[1], sorted_pt_inds[i]);
        }
    }

    i = 2;
    {
        Vector3d p0{point(points, sorted_pt_inds[0])};
        const auto v1 = (Vector3d{point(points, sorted_pt_inds[1])} - p0).eval();
        while (i < n_points) {
            Vector3d pi{point(points, sorted_pt_inds[i])};
            const auto v2 = pi - p0;
            const auto n = v1.cross(v2);
            if (n.norm() / bbox_size2 < epsilon) {
                i += 1;
            } else {
                break;
            }
        }
        if (i > 2) {
            if (i == n_points) return {};
            std::swap(sorted_pt_inds[2], sorted_pt_inds[i]);
        }
    }
    i = 3;
    double ori = 0.0;
    while (i < n_points) {
        ori = orient3dfast(
            point(points, sorted_pt_inds[0]), point(points, sorted_pt_inds[1]), point(points, sorted_pt_inds[2]),
            point(points, sorted_pt_inds[i])
        );
        if (std::fabs(ori) / bbox_size3 < epsilon) {
            i += 1;
        } else {
            break;
        }
    }
    if (i > 3) {
        if (i == n_points) return {};
        std::swap(sorted_pt_inds[3], i);
    }

    // follow the right rule
    if (ori > 0.0) {
        std::swap(sorted_pt_inds[0], sorted_pt_inds[1]);
    }

    Tetrahedrons tets{points, n_points, {}, std::vector<uint32_t>(n_points + 1, INVALID)};
    TriFace search_tet =
        initial_delaunay(tets, sorted_pt_inds[0], sorted_pt_inds[1], sorted_pt_inds[2], sorted_pt_inds[3]);
    for (i = 4; i < n_points; i++) {
        if (!insert_vertex_bw(tets, sorted_pt_inds[i], search_tet)) {
            break;
        }
    }
    return tets;
}
