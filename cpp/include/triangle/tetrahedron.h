#pragma once

#include <predicates/predicates.h>

#include <array>
#include <limits>
#include <memory>
#include <vector>

static constexpr std::array<uint32_t, 12> enext_table() noexcept {
    std::array<uint32_t, 12> enexttbl;
    for (uint32_t i = 0; i < 12; i++) {
        enexttbl[i] = (i + 4) % 12;
    }
    return enexttbl;
}

static constexpr std::array<uint32_t, 12> eprev_table() noexcept {
    std::array<uint32_t, 12> eprevtbl;
    for (uint32_t i = 0; i < 12; i++) {
        eprevtbl[i] = (i + 8) % 12;
    }
    return eprevtbl;
}

static constexpr std::array<uint32_t, 12>
enext_esym_table(const std::array<uint32_t, 12>& esymtbl, const std::array<uint32_t, 12>& enexttbl) noexcept {
    std::array<uint32_t, 12> enextesymtbl;
    for (uint32_t i = 0; i < 12; i++) {
        enextesymtbl[i] = esymtbl[enexttbl[i]];
    }
    return enextesymtbl;
}

static constexpr std::array<uint32_t, 12>
eprev_esym_table(const std::array<uint32_t, 12>& esymtbl, const std::array<uint32_t, 12>& eprevtbl) noexcept {
    std::array<uint32_t, 12> eprevesymtbl;
    for (uint32_t i = 0; i < 12; i++) {
        eprevesymtbl[i] = esymtbl[eprevtbl[i]];
    }
    return eprevesymtbl;
}

static constexpr std::array<std::array<uint32_t, 12>, 12> bond_table() noexcept {
    std::array<std::array<uint32_t, 12>, 12> bondtbl;
    for (uint32_t i = 0; i < 12; i++) {
        for (uint32_t j = 0; j < 12; j++) {
            bondtbl[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
        }
    }
    return bondtbl;
}

static constexpr std::array<std::array<uint32_t, 12>, 12> fsym_table() noexcept {
    std::array<std::array<uint32_t, 12>, 12> fsymtbl;
    for (uint32_t i = 0; i < 12; i++) {
        for (uint32_t j = 0; j < 12; j++) {
            fsymtbl[i][j] = (j + 12 - (i & 12)) % 12;
        }
    }
    return fsymtbl;
}

static constexpr std::array<uint32_t, 12> face_pivot1(const std::array<uint32_t, 12>& esymtbl) noexcept {
    std::array<uint32_t, 12> facepivot1;
    for (uint32_t i = 0; i < 12; i++) {
        facepivot1[i] = (esymtbl[i] & 3);
    }
    return facepivot1;
}

static constexpr std::array<std::array<uint32_t, 12>, 12>
face_pivot2(const std::array<std::array<uint32_t, 12>, 12>& fsymtbl, const std::array<uint32_t, 12>& esymtbl) noexcept {
    std::array<std::array<uint32_t, 12>, 12> facepivot2;
    for (uint32_t i = 0; i < 12; i++) {
        for (uint32_t j = 0; j < 12; j++) {
            facepivot2[i][j] = fsymtbl[esymtbl[i]][j];
        }
    }
    return facepivot2;
}

struct TriFace {
    static constexpr std::array<uint32_t, 12> enexttbl{enext_table()};
    static constexpr std::array<uint32_t, 12> eprevtbl{eprev_table()};
    static constexpr std::array<uint32_t, 12> esymtbl{{9, 6, 11, 4, 3, 7, 1, 5, 10, 0, 8, 2}};
    static constexpr std::array<uint32_t, 12> enextesymtbl{enext_esym_table(esymtbl, enexttbl)};
    static constexpr std::array<uint32_t, 12> eprevesymtbl{eprev_esym_table(esymtbl, eprevtbl)};
    static constexpr uint32_t INVALID = std::numeric_limits<uint32_t>::max();

    uint32_t tet{INVALID};
    uint32_t ver{11}; // version
    TriFace() {}
    TriFace(uint32_t&& t, uint32_t&& v) : tet{t}, ver(v) {}
    TriFace(const uint32_t& t, const uint32_t& v) : tet(t), ver(v) {}
    TriFace(const TriFace& f) : tet(f.tet), ver(f.ver) {}

    TriFace& operator=(const TriFace& f) {
        tet = f.tet;
        ver = f.ver;
        return *this;
    }
    void set(const uint32_t t, const uint32_t v) {
        tet = t;
        ver = v;
    }

    // next edge
    static void enext(const TriFace& t1, TriFace& t2) {
        t2.tet = t1.tet;
        t2.ver = enexttbl[t1.ver];
    }
    void enext_self() { ver = enexttbl[ver]; }

    // prev edge
    static void eprev(const TriFace& t1, TriFace& t2) {
        t2.tet = t1.tet;
        t2.ver = eprevtbl[t1.ver];
    }
    void eprev_self() { ver = eprevtbl[ver]; }

    // symmetric edge
    static void esym(const TriFace& t1, TriFace& t2) {
        t2.tet = t1.tet;
        t2.ver = esymtbl[t1.ver];
    }
    void esym_self() { ver = esymtbl[ver]; }

    // next sym edge
    static void enext_esym(const TriFace& t1, TriFace& t2) {
        t2.tet = t1.tet;
        t2.ver = enextesymtbl[t1.ver];
    }
    void enext_esym_self() { ver = enextesymtbl[ver]; }

    // prev sym edge
    static void eprev_esym(const TriFace& t1, TriFace& t2) {
        t2.tet = t1.tet;
        t2.ver = eprevesymtbl[t1.ver];
    }
    void eprev_esym_self() { ver = eprevesymtbl[ver]; }
};

struct Tet {
    std::array<uint32_t, 4> data; // tet vertices
    std::array<TriFace, 4> nei;   // neighbor
    uint8_t mask{0};

    Tet(std::array<uint32_t, 4>&& d, std::array<TriFace, 4>&& n) : data(std::move(d)), nei(std::move(n)), mask{0} {}

    uint32_t index(const uint32_t vid) const {
        uint32_t i = 0;
        for (; i < 4; i++) {
            if (data[i] == vid) {
                return i;
            }
        }
        return i;
    }
};

struct TetMesh {
    static constexpr std::array<std::array<uint32_t, 12>, 12> bondtbl{bond_table()};
    static constexpr std::array<std::array<uint32_t, 12>, 12> fsymtbl{fsym_table()};
    static constexpr std::array<uint32_t, 12> facepivot1{face_pivot1(TriFace::esymtbl)};
    static constexpr std::array<std::array<uint32_t, 12>, 12> facepivot2{face_pivot2(fsymtbl, TriFace::esymtbl)};
    static constexpr std::array<uint32_t, 12> orgpivot{{3, 3, 1, 1, 2, 0, 0, 2, 1, 2, 3, 0}};
    static constexpr std::array<uint32_t, 12> destpivot{{2, 0, 0, 2, 1, 2, 3, 0, 3, 3, 1, 1}};
    static constexpr std::array<uint32_t, 12> apexpivot{{1, 2, 3, 0, 3, 3, 1, 1, 2, 0, 0, 2}};
    static constexpr std::array<uint32_t, 12> oppopivot{{0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3}};
    static constexpr std::array<uint32_t, 4> vpivot{{11, 8, 9, 10}};

    const double* points{nullptr};
    const uint32_t n_points{0};
    std::vector<Tet> tets;
    std::vector<uint32_t> p2t; // point to tet
    static TetMesh tetrahedralize(const double* points, const uint32_t n_points, const double epsilon);

    bool is_hull_tet(const uint32_t idx) const { return tets[idx].data[3] == n_points; }

    void bond(const TriFace& t1, const TriFace& t2) {
        TriFace& f1 = tets[t1.tet].nei[t1.ver & 3];
        f1.tet = t2.tet;
        f1.ver = bondtbl[t1.ver][t2.ver];

        TriFace& f2 = tets[t2.tet].nei[t2.ver & 3];
        f2.tet = t1.tet;
        f2.ver = bondtbl[t2.ver][t1.ver];
    }
    void fsym(const TriFace& t1, TriFace& t2) const {
        const TriFace& nf = tets[t1.tet].nei[t1.ver & 3];
        t2.tet = nf.tet;
        t2.ver = fsymtbl[t1.ver][nf.ver];
    }

    void fsym_self(TriFace& t) const {
        const TriFace& nf = tets[t.tet].nei[t.ver & 3];
        t.tet = nf.tet;
        t.ver = fsymtbl[t.ver][nf.ver];
    }

    void fnext_self(TriFace& t) const {
        const TriFace& nf = tets[t.tet].nei[facepivot1[t.ver]];
        t.tet = nf.tet;
        t.ver = facepivot2[t.ver][nf.ver];
    }
    void infect(const uint32_t t) { tets[t].mask |= 1; }

    void uninfect(const uint32_t t) { tets[t].mask &= ~1; }

    bool infected(const uint32_t t) const { return (tets[t].mask & 1) != 0; }

    void mark_test(const uint32_t t) { tets[t].mask |= 2; }

    void unmark_test(const uint32_t t) { tets[t].mask &= ~2; }

    bool mark_tested(const uint32_t t) const { return (tets[t].mask & 2) != 0; }

    uint32_t org(const TriFace& f) const { return tets[f.tet].data[orgpivot[f.ver]]; }

    uint32_t dest(const TriFace& f) const { return tets[f.tet].data[destpivot[f.ver]]; }

    uint32_t apex(const TriFace& f) const { return tets[f.tet].data[apexpivot[f.ver]]; }

    uint32_t oppo(const TriFace& f) const { return tets[f.tet].data[oppopivot[f.ver]]; }

    const double* point(const uint32_t idx) const { return &points[idx * 3]; }

    // tets incident to vertex
    uint32_t incident(const uint32_t vid, std::vector<uint32_t>& result) {
        uint32_t count = 0;
        const auto push = [&count, &result, this](const uint32_t t) {
            if (count < result.size()) {
                result[count] = t;
            } else {
                result.emplace_back(t);
            }
            mark_test(t);
            count += 1;
        };
        push(p2t[vid]);
        for (uint32_t i = 0; i < count; i++) {
            const uint32_t t = result[i];
            const uint32_t ind = tets[t].index(vid);
            for (uint32_t j = 0; j < 4; j++) {
                if (j == ind) {
                    continue;
                }
                const uint32_t nei = tets[t].nei[j].tet;
                if (!mark_tested(nei) && !is_hull_tet(nei)) {
                    push(nei);
                }
            }
        }
        for (uint32_t i = 0; i < count; i++) {
            unmark_test(result[i]);
        }
        return count;
    }
    // tets incident to edge [va, vb]
    uint32_t incident(const TriFace& edge, std::vector<uint32_t>& result) {
        TriFace spin{edge};
        uint32_t count = 0;
        do {
            if (!is_hull_tet(spin.tet)) {
                if (count < result.size()) {
                    result[count++] = spin.tet;
                } else {
                    count += 1;
                    result.emplace_back(spin.tet);
                }
            }
            fnext_self(spin);
        } while (spin.tet != edge.tet);
        return count;
    }

    int orient3d(const uint32_t pa, const uint32_t pb, const uint32_t pc, const uint32_t pd) const {
        const double ret = ::orient3d(&points[pa * 3], &points[pb * 3], &points[pc * 3], &points[pd * 3]);
        return (ret > 0.0) - (ret < 0.0);
    }
};
