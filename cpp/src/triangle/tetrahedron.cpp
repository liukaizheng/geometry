#include <predicates/predicates.h>
#include <triangle/tetrahedron.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <unordered_map>
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

constexpr std::array<uint32_t, 12> epivot{{4, 5, 2, 11, 4, 5, 2, 11, 4, 5, 2, 11}};

// Insphere test with symbolic perturbation
inline int insphere_s(
    const TetMesh& tets, const uint32_t pa, const uint32_t pb, const uint32_t pc, const uint32_t pd, const uint32_t pe
) {
    const double* points = tets.points;
    const double sign = insphere(&points[pa * 3], &points[pb * 3], &points[pc * 3], &points[pd * 3], &points[pe * 3]);
    if (sign != 0.0) {
        return (sign > 0.0) - (sign < 0.0);
    }

    uint32_t pt[5] = {pa, pb, pc, pd, pe};

    int swaps = 0; // Record the total number of swaps.
    int n = 5;
    int count = 0;
    do {
        count = 0;
        n = n - 1;
        for (int i = 0; i < n; i++) {
            if (pt[i] > pt[i + 1]) {
                std::swap(pt[i], pt[i + 1]);
                count++;
            }
        }
        swaps += count;
    } while (count > 0); // Continue if some points are swapped.

    int oriA = tets.orient3d(pt[1], pt[2], pt[3], pt[4]);
    if (oriA != 0.0) {
        // Flip the sign if there are odd number of swaps.
        if ((swaps % 2) != 0) oriA = -oriA;
        return oriA;
    }

    int oriB = -tets.orient3d(pt[0], pt[2], pt[3], pt[4]);
    // Flip the sign if there are odd number of swaps.
    if ((swaps % 2) != 0) oriB = -oriB;
    return oriB;
}

inline void
make_tet(TetMesh& tets, TriFace& face, const uint32_t pa, const uint32_t pb, const uint32_t pc, const uint32_t pd) {
    face.tet = static_cast<uint32_t>(tets.tets.size());
    face.ver = 11;
    tets.tets.emplace_back(std::array<uint32_t, 4>{{pa, pb, pc, pd}}, std::array<TriFace, 4>{{{}, {}, {}, {}}});
}

inline TriFace
initial_delaunay(TetMesh& tets, const uint32_t pa, const uint32_t pb, const uint32_t pc, const uint32_t pd) {
    TriFace firsttet, tetopa, tetopb, tetopc, tetopd;

    make_tet(tets, firsttet, pa, pb, pc, pd);
    make_tet(tets, tetopa, pb, pc, pd, tets.n_points);
    make_tet(tets, tetopb, pc, pa, pd, tets.n_points);
    make_tet(tets, tetopc, pa, pb, pd, tets.n_points);
    make_tet(tets, tetopd, pb, pa, pc, tets.n_points);

    // Connect hull tetrahedra to firsttet (at four faces of firsttet).
    TriFace worktet;
    tets.bond(firsttet, tetopd);
    TriFace::esym(firsttet, worktet);
    tets.bond(worktet, tetopc); // ab
    TriFace::enext_esym(firsttet, worktet);
    tets.bond(worktet, tetopa); // bc
    TriFace::eprev_esym(firsttet, worktet);
    tets.bond(worktet, tetopb); // ca

    // Connect hull tetrahedra together (at six edges of firsttet).
    TriFace worktet1;
    TriFace::esym(tetopc, worktet);
    TriFace::esym(tetopd, worktet1);
    tets.bond(worktet, worktet1); // ab
    TriFace::esym(tetopa, worktet);
    TriFace::eprev_esym(tetopd, worktet1);
    tets.bond(worktet, worktet1); // bc
    TriFace::esym(tetopb, worktet);
    TriFace::enext_esym(tetopd, worktet1);
    tets.bond(worktet, worktet1); // ca
    TriFace::eprev_esym(tetopc, worktet);
    TriFace::enext_esym(tetopb, worktet1);
    tets.bond(worktet, worktet1); // da
    TriFace::eprev_esym(tetopa, worktet);
    TriFace::enext_esym(tetopc, worktet1);
    tets.bond(worktet, worktet1); // db
    TriFace::eprev_esym(tetopb, worktet);
    TriFace::enext_esym(tetopa, worktet1);
    tets.bond(worktet, worktet1); // dc

    tets.p2t[pa] = 0;
    tets.p2t[pb] = 0;
    tets.p2t[pc] = 0;
    tets.p2t[pd] = 0;

    tets.p2t[tets.n_points] = tetopa.tet;

    return firsttet;
}

enum class LocateResult { OUTSIDE, ONVERTEX, ONEDGE, ONFACE, INTETRAHEDRON };

inline LocateResult locate_dt(const TetMesh& tets, const uint32_t pid, TriFace& searchtet) {
    if (tets.is_hull_tet(searchtet.tet)) {
        searchtet.tet = tets.tets[searchtet.tet].nei[3].tet;
    }

    for (searchtet.ver = 0; searchtet.ver < 4; searchtet.ver++) {
        const int ori = tets.orient3d(tets.org(searchtet), tets.dest(searchtet), tets.apex(searchtet), pid);
        if (ori < 0) break;
    }

    LocateResult loc = LocateResult::OUTSIDE;
    while (true) {
        uint32_t toppo = tets.oppo(searchtet);

        // Check if the vertex is we seek.
        if (toppo == pid) {
            // Adjust the origin of searchtet to be searchpt.
            searchtet.esym_self();
            searchtet.eprev_self();
            loc = LocateResult::ONVERTEX; // return ONVERTEX;
            break;
        }

        const int oriorg = tets.orient3d(tets.dest(searchtet), tets.apex(searchtet), toppo, pid);
        if (oriorg < 0) {
            searchtet.enext_esym_self();
        } else {
            const int oridest = tets.orient3d(tets.apex(searchtet), tets.org(searchtet), toppo, pid);
            if (oridest < 0) {
                searchtet.eprev_esym_self();
            } else {
                const int oriapex = tets.orient3d(tets.org(searchtet), tets.dest(searchtet), toppo, pid);
                if (oriapex < 0) {
                    searchtet.esym_self();
                } else {
                    // oriorg >= 0, oridest >= 0, oriapex >= 0 ==> found the point.
                    // The point we seek must be on the boundary of or inside this
                    //   tetrahedron. Check for boundary cases first.
                    if (oriorg == 0) {
                        // Go to the face opposite to origin.
                        searchtet.enext_esym_self();
                        if (oridest == 0) {
                            searchtet.eprev_self(); // edge oppo->apex
                            if (oriapex == 0) {
                                // oppo is duplicated with p.
                                loc = LocateResult::ONVERTEX; // return ONVERTEX;
                                break;
                            }
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        if (oriapex == 0) {
                            searchtet.enext_self();     // edge dest->oppo
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        loc = LocateResult::ONFACE; // return ONFACE;
                        break;
                    }
                    if (oridest == 0) {
                        // Go to the face opposite to destination.
                        searchtet.eprev_esym_self();
                        if (oriapex == 0) {
                            searchtet.eprev_self();     // edge oppo->org
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        loc = LocateResult::ONFACE; // return ONFACE;
                        break;
                    }
                    if (oriapex == 0) {
                        // Go to the face opposite to apex
                        searchtet.esym_self();
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
inline bool insert_vertex_bw(TetMesh& tets, const uint32_t pid, TriFace& searchtet) {
    const LocateResult loc = locate_dt(tets, pid, searchtet);
    std::vector<uint32_t> cave_oldtet_list;
    if (loc == LocateResult::OUTSIDE) {
        tets.infect(searchtet.tet);
        cave_oldtet_list.emplace_back(searchtet.tet);
    } else if (loc == LocateResult::INTETRAHEDRON) {
        tets.infect(searchtet.tet);
        cave_oldtet_list.emplace_back(searchtet.tet);
    } else if (loc == LocateResult::ONFACE) {
        tets.infect(searchtet.tet);
        cave_oldtet_list.emplace_back(searchtet.tet);
        const uint32_t nei_tet = tets.tets[searchtet.tet].nei[searchtet.ver & 3].tet;
        tets.infect(nei_tet);
        cave_oldtet_list.emplace_back(nei_tet);
    } else if (loc == LocateResult::ONEDGE) {
        TriFace spintet = searchtet;
        while (true) {
            tets.infect(spintet.tet);
            cave_oldtet_list.emplace_back(spintet.tet);
            tets.fnext_self(spintet);
            if (spintet.tet == searchtet.tet) break;
        }
    } else if (loc == LocateResult::ONVERTEX) {
        // The point already exist. Do nothing and return.
        return false;
    }
    uint32_t cavetid, neightid;
    std::vector<TriFace> cave_bdry_list;
    for (uint32_t i = 0; i < cave_oldtet_list.size(); i++) {
        cavetid = cave_oldtet_list[i];
        for (uint32_t ver = 0; ver < 4; ver++) {
            neightid = tets.tets[cavetid].nei[ver].tet;
            if (tets.infected(neightid)) continue;
            bool enqflag = false;
            if (!tets.mark_tested(neightid)) {
                const auto& pts = tets.tets[neightid].data;
                if (!tets.is_hull_tet(neightid)) {
                    enqflag = insphere_s(tets, pts[0], pts[1], pts[2], pts[3], pid) < 0;
                } else {
                    const int ori = tets.orient3d(pts[0], pts[1], pts[2], pid);
                    if (ori < 0) {
                        enqflag = true;
                    } else if (ori == 0) {
                        uint32_t neineitet = tets.tets[neightid].nei[3].tet;
                        const auto& nei_pts = tets.tets[neineitet].data;
                        enqflag = insphere_s(tets, nei_pts[0], nei_pts[1], nei_pts[2], nei_pts[3], pid) < 0;
                    }
                }
                tets.mark_test(neightid);
            }

            if (enqflag) {
                tets.infect(neightid);
                cave_oldtet_list.emplace_back(neightid);
            } else {
                // A boundary face.
                cave_bdry_list.emplace_back(cavetid, ver);
            }
        }
    }

    constexpr uint32_t row_v08_tbl[12] = {8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7};
    constexpr uint32_t row_v11_tbl[12] = {8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7};
    constexpr uint32_t col_v01_tbl[12] = {1, 1, 1, 1, 5, 5, 5, 5, 9, 9, 9, 9};
    constexpr uint32_t col_v02_tbl[12] = {2, 2, 2, 2, 6, 6, 6, 6, 10, 10, 10, 10};
    constexpr uint32_t col_v08_tbl[12] = {8, 8, 8, 8, 0, 0, 0, 0, 4, 4, 4, 4};
    constexpr uint32_t col_v11_tbl[12] = {11, 11, 11, 11, 3, 3, 3, 3, 7, 7, 7, 7};
        
    const uint32_t f_out = static_cast<uint32_t>(cave_bdry_list.size());
    const uint32_t v_out = (f_out + 4) >> 1;
        
    static std::array<TriFace, 4096> bw_faces;
    TriFace* tmp_bw_faces = nullptr;
    uint32_t shiftbits = 0;
    std::unique_ptr<TriFace[]> dyn_bw_faces;
    if (v_out < 64) {
        shiftbits = 6;
        tmp_bw_faces = bw_faces.data();
    } else if (v_out < 1024) {
        // Dynamically allocate an array to store the adjacencies.
        uint32_t tmp = v_out;
        shiftbits = 1;
        while ((tmp >>= 1)) shiftbits++;
        const uint32_t arysize = 1 << shiftbits;
        dyn_bw_faces = std::make_unique<TriFace[]>(arysize * arysize);
        tmp_bw_faces = dyn_bw_faces.get();
    }

    uint32_t local_vcount = 0;
    if (v_out < 1024) {
        // pid to local vertex id
        std::unordered_map<uint32_t, uint32_t> pmap;
        for (uint32_t i = 0; i < f_out; i++) {
            TriFace& oldtet = cave_bdry_list[i];
            // Get the tet outside the cavity.
            TriFace neightet = tets.tets[oldtet.tet].nei[oldtet.ver];
            tets.unmark_test(neightet.tet);

            if (tets.is_hull_tet(oldtet.tet)) {
                // neightet.tet may be also a hull tet (=> oldtet is a hull edge).
                neightet.ver = epivot[neightet.ver];
            }

            // Create a new tet in the cavity.
            uint32_t v[3] = {tets.dest(neightet), tets.org(neightet), tets.apex(neightet)};
            TriFace newtet;
            make_tet(tets, newtet, v[1], v[0], pid, v[2]);
            tets.tets[newtet.tet].nei[2] = neightet;
            tets.tets[neightet.tet].nei[neightet.ver & 3].set(newtet.tet, col_v02_tbl[neightet.ver]);

            uint32_t sidx[3];
            // Fill the adjacency matrix, and count v_out.
            for (uint32_t j = 0; j < 3; j++) {
                const uint32_t tid = tets.p2t[v[j]];
                if (tets.tets[tid].data[2] != pid) {
                    pmap.emplace(v[j], local_vcount++);
                    tets.p2t[v[j]] = newtet.tet;
                }
                sidx[j] = pmap[v[j]];
            }

            neightet.tet = newtet.tet;
            neightet.ver = 11;
            tmp_bw_faces[(sidx[1] << shiftbits) | sidx[0]] = neightet;
            neightet.ver = 1;
            tmp_bw_faces[(sidx[2] << shiftbits) | sidx[1]] = neightet;
            neightet.ver = 8;
            tmp_bw_faces[(sidx[0] << shiftbits) | sidx[2]] = neightet;

            oldtet = newtet;
        }

        // randomly pick a new tet
        searchtet = cave_bdry_list[f_out >> 1];
        tets.p2t[pid] = searchtet.tet;

        for (uint32_t i = 0; i < f_out; i++) {
            TriFace neightet = cave_bdry_list[i];
            if (tets.tets[neightet.tet].nei[3].tet == INVALID) {
                neightet.ver = 11;
                const uint32_t j = pmap[tets.org(neightet)];
                const uint32_t k = pmap[tets.dest(neightet)];
                const TriFace& neineitet = tmp_bw_faces[(k << shiftbits) | j];
                tets.tets[neightet.tet].nei[3].set(neineitet.tet, row_v11_tbl[neineitet.ver]);
                tets.tets[neineitet.tet].nei[neineitet.ver & 3].set(neightet.tet, col_v11_tbl[neineitet.ver]);
            }
            if (tets.tets[neightet.tet].nei[1].tet == INVALID) {
                neightet.ver = 1;
                const uint32_t j = pmap[tets.org(neightet)];
                const uint32_t k = pmap[tets.dest(neightet)];
                const TriFace& neineitet = tmp_bw_faces[(k << shiftbits) | j];
                tets.tets[neightet.tet].nei[1].set(neineitet.tet, neineitet.ver);
                tets.tets[neineitet.tet].nei[neineitet.ver & 3].set(neightet.tet, col_v01_tbl[neineitet.ver]);
            }
            if (tets.tets[neightet.tet].nei[0].tet == INVALID) {
                neightet.ver = 8;
                const uint32_t j = pmap[tets.org(neightet)];
                const uint32_t k = pmap[tets.dest(neightet)];
                const TriFace& neineitet = tmp_bw_faces[(k << shiftbits) | j];
                tets.tets[neightet.tet].nei[0].set(neineitet.tet, row_v08_tbl[neineitet.ver]);
                tets.tets[neineitet.tet].nei[neineitet.ver & 3].set(neightet.tet, col_v08_tbl[neineitet.ver]);
            }
        }
    } else {
        // Fill a very large cavity with original neighboring searching method.
        for (uint32_t i = 0; i < f_out; i++) {
            TriFace& oldtet = cave_bdry_list[i];
            TriFace neightet = tets.tets[oldtet.tet].nei[oldtet.ver];

            tets.unmark_test(neightet.tet);
            if (tets.is_hull_tet(oldtet.tet)) {
                // neightet.tet may be also a hull tet (=> oldtet is a hull edge).
                neightet.ver = epivot[neightet.ver];
            }

            uint32_t v[3] = {tets.dest(neightet), tets.org(neightet), tets.apex(neightet)};
            TriFace newtet;
            make_tet(tets, newtet, v[1], v[0], pid, v[2]);

            tets.tets[newtet.tet].nei[2].set(neightet.tet, neightet.ver);
            tets.tets[neightet.tet].nei[neightet.ver & 3].set(newtet.tet, col_v02_tbl[neightet.ver]);

            // Fill the adjacency matrix, and count v_out.
            for (uint32_t j = 0; j < 3; j++) {
                const uint32_t tid = tets.p2t[v[j]];
                if (tets.tets[tid].data[2] != pid) {
                    local_vcount += 1;
                    tets.p2t[v[j]] = newtet.tet;
                }
            }
        }

        // randomly pick a new tet
        searchtet = cave_bdry_list[f_out >> 1];
        tets.p2t[pid] = searchtet.tet;

        for (uint32_t i = 0; i < f_out; i++) {
            TriFace oldtet = cave_bdry_list[i];

            TriFace neightet, newtet;
            tets.fsym(oldtet, neightet);
            tets.fsym(neightet, newtet);
            // Oldtet and newtet must be at the same directed edge.
            // Connect the three other faces of this newtet.
            for (uint32_t j = 0; j < 3; j++) {
                TriFace::esym(newtet, neightet); // Go to the face.
                if (tets.tets[neightet.tet].nei[neightet.ver & 3].tet == INVALID) {
                    // Find the adjacent face of this newtet.
                    TriFace spintet = oldtet;
                    while (true) {
                        tets.fnext_self(spintet);
                        if (!tets.infected(spintet.tet)) break;
                    }
                    TriFace neineitet;
                    tets.fsym(spintet, neineitet);
                    neineitet.esym_self();
                    tets.bond(neightet, neineitet);
                }
                newtet.enext_self();
                oldtet.enext_self();
            }
        }
    }

    for (uint32_t i = 0; i < cave_oldtet_list.size(); i++) {
        tets.tets[cave_oldtet_list[i]].mask = static_cast<uint8_t>(-1);
    }
    return true;
}

TetMesh TetMesh::tetrahedralize(const double* points, const uint32_t n_points, const double epsilon) {
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
        std::swap(sorted_pt_inds[3], sorted_pt_inds[i]);
    }

    // follow the right rule
    if (ori > 0.0) {
        std::swap(sorted_pt_inds[0], sorted_pt_inds[1]);
    }

    initstaticfilter(
        std::fmax(std::fabs(min_corner[0]), std::fabs(max_corner[0])),
        std::fmax(std::fabs(min_corner[1]), std::fabs(max_corner[1])),
        std::fmax(std::fabs(min_corner[2]), std::fabs(max_corner[2]))
    );

    TetMesh tets{points, n_points, {}, std::vector<uint32_t>(n_points + 1, INVALID)};
    TriFace search_tet =
        initial_delaunay(tets, sorted_pt_inds[0], sorted_pt_inds[1], sorted_pt_inds[2], sorted_pt_inds[3]);
    for (i = 4; i < n_points; i++) {
        if (!insert_vertex_bw(tets, sorted_pt_inds[i], search_tet)) {
            break;
        }
    }
    auto it = std::remove_if(tets.tets.begin(), tets.tets.end(), [](const Tet& t) {
        return t.mask == static_cast<uint8_t>(-1);
    });
    tets.tets.erase(it, tets.tets.end());
    return tets;
}
