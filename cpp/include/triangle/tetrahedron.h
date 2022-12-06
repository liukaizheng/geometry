#pragma once

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

struct TriFace {
    static constexpr std::array<uint32_t, 12> enexttbl{enext_table()};
    static constexpr std::array<uint32_t, 12> eprevtbl{eprev_table()};
    static constexpr std::array<uint32_t, 12> esymtbl{{9, 6, 11, 4, 3, 7, 1, 5, 10, 0, 8, 2}};
    static constexpr std::array<uint32_t, 12> enextesymtbl{enext_esym_table(esymtbl, enexttbl)};
    static constexpr std::array<uint32_t, 12> eprevesymtbl{eprev_esym_table(esymtbl, eprevtbl)};

    uint32_t tet{std::numeric_limits<uint32_t>::max()};
    uint32_t ver{11}; // version
    TriFace() {}
    TriFace(uint32_t&& t, uint32_t&& v) : tet{t}, ver(v) {}
    TriFace(const uint32_t& t, const uint32_t& v) : tet(t), ver(v) {}

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
};

static constexpr std::array<std::array<uint32_t, 12>, 12> bond_table() noexcept {
    std::array<std::array<uint32_t, 12>, 12> bondtbl;
    for (uint32_t i = 0; i < 12; i++) {
        for (uint32_t j = 0; j < 12; j++) {
            bondtbl[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
        }
    }
    return bondtbl;
}

struct TetMesh {
    static constexpr std::array<std::array<uint32_t, 12>, 12> bondtbl{bond_table()};
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
};
