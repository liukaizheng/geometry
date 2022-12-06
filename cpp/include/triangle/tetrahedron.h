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

struct TetMesh {
    const double* points { nullptr };
    const uint32_t n_points { 0 };
    std::vector<Tet> tets;
    std::vector<uint32_t> p2t;    // point to tet
    static TetMesh tetrahedralize(const double* points, const uint32_t n_points, const double epsilon);
    
    bool is_hull_tet(const uint32_t idx) const { return tets[idx].data[3] == n_points; }
};
