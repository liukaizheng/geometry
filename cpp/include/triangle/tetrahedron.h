#pragma once

#include <array>
#include <vector>
#include <limits>
#include <memory>

struct TriFace {
    uint32_t tet { std::numeric_limits<uint32_t>::max() };
    uint32_t ver { 11 }; // version
    TriFace() {}
    TriFace(uint32_t&& t, uint32_t && v): tet{t}, ver(v) {}
    TriFace(const uint32_t& t, const uint32_t& v): tet(t), ver(v) {}
    
    TriFace& operator=(const TriFace& f) {
        tet = f.tet;
        ver = f.ver;
        return *this;
    }
    void set(const uint32_t t, const uint32_t v) {
        tet = t;
        ver = v;
    }
};

struct Tet {
    std::array<uint32_t, 4> data; // tet vertices
    std::array<TriFace, 4> nei;   // neighbor  
    uint8_t mask {0};
    Tet(std::array<uint32_t, 4>&& d, std::array<TriFace, 4>&& n): data(std::move(d)), nei(std::move(n)), mask{0} {}
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
    const double* points { nullptr };
    const uint32_t n_points { 0 };
    std::vector<Tet> tets;
    std::vector<uint32_t> p2t;    // point to tet
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
