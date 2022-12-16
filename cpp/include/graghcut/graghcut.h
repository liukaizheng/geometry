#pragma once

#include <vector>

class GraphCut {
    static constexpr uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    struct Arc {
        uint32_t head;          // head
        uint32_t next{INVALID}; // next
        uint32_t sister;        // sister
        double r_cap{0};        // residual capacity
        Arc() {}
        Arc(const uint32_t h, const uint32_t n, const uint32_t s, const double c)
            : head(h), next(n), sister(s), r_cap(c) {}
    };

    const uint32_t n_nodes;
    const uint32_t TERMINAL;
    const uint32_t ORPHAN;
    uint32_t time{0};
    uint32_t queue_first[2]{INVALID, INVALID};
    uint32_t queue_last[2]{INVALID, INVALID};
    uint32_t orphan_first{INVALID};
    uint32_t orphan_last{INVALID};
    std::vector<double> tr_cap; // residual capacities of nodes
    double flow{0.0};
    std::vector<Arc> arcs;             // all arcs
    std::vector<uint32_t> dist;        // distance to the terminal
    std::vector<uint32_t> first_arc;   // first arcs of nodes
    std::vector<uint32_t> next;        // next active node
    std::vector<uint32_t> next_orphan; // next orphan node
    std::vector<uint32_t> parent;      // parent arcs of nodes
    std::vector<uint32_t> ts;          // time stamp showing when dist was computed
    std::vector<bool> is_sink;         // whether the node is in the source

  public:
    GraphCut(const uint32_t n, double* source_cap, double* sink_cap);

    void set_active(const uint32_t i) {
        if (queue_last[1] == INVALID) {
            queue_first[1] = i;
        } else {
            next[queue_last[1]] = i;
        }
        queue_last[1] = i;
        next[i] = i;
    }

    uint32_t next_active() {
        while (true) {
            uint32_t i = queue_first[0];
            if (i == INVALID) {
                queue_first[0] = i = queue_first[1];
                queue_last[0] = queue_last[1];
                queue_first[1] = INVALID;
                queue_last[1] = INVALID;
                if (i == INVALID) {
                    return INVALID;
                }
            }
            if (next[i] == i) {
                queue_first[0] = queue_last[0] = INVALID;
            } else {
                queue_first[0] = next[i];
            }
            next[i] = INVALID;
            if (parent[i] != INVALID) {
                return i;
            }
        }
    }

    void add_edge(const uint32_t i, const uint32_t j, const double cap, const double rev_cap) {
        const uint32_t arc_id = static_cast<uint32_t>(arcs.size());
        const uint32_t sisiter_arc_id = arc_id + 1;
        arcs.emplace_back(j, first_arc[i], sisiter_arc_id, cap); // arc
        arcs.emplace_back(i, first_arc[j], arc_id, rev_cap);     // sisiter arc
    }
    void augment(Arc* minnle_arc);
    double max_flow();
};
