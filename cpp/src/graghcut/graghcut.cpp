#include <graghcut/graghcut.h>

GraphCut::GraphCut(const uint32_t n, double* source_cap, double* sink_cap)
    : n_nodes(n), TERMINAL(n * 2), ORPHAN(TERMINAL + 1) {
    dist.resize(n, 1);
    is_sink.resize(n);
    next.resize(n, INVALID);
    next_orphan.resize(n, INVALID);
    parent.resize(n, TERMINAL);
    ts.resize(n, 0);
    for (uint32_t i = 0; i < n; i++) {
        tr_cap[i] = source_cap[i] - sink_cap[i];
        flow += source_cap[i] < sink_cap[i] ? source_cap[i] : sink_cap[i];
        is_sink[i] = tr_cap[i] < 0.0;
        if (tr_cap[i] == 0.0) {
            parent[i] = INVALID;
        } else {
            set_active(i);
        }
    }
    first_arc.resize(n, INVALID);
}

void GraphCut::augment(Arc* middle_arc) {
    double bottle_neck = middle_arc->r_cap;
    uint32_t i = arcs[middle_arc->sister].head;
    uint32_t aid = INVALID;
    while (true) {
        aid = parent[i];
        if (aid == TERMINAL) {
            break;
        }
        const Arc& sister = arcs[arcs[aid].sister];
        if (bottle_neck > sister.r_cap) {
            bottle_neck = sister.r_cap;
        }
    }

    if (bottle_neck > tr_cap[i]) {
        bottle_neck = tr_cap[i];
    }

    i = middle_arc->head;
    while (true) {
        aid = parent[i];
        if (aid == TERMINAL) {
            break;
        }
        const Arc& arc = arcs[aid];
        if (bottle_neck > arc.r_cap) {
            bottle_neck = arc.r_cap;
        }
    }

    if (bottle_neck > -tr_cap[i]) {
        bottle_neck = -tr_cap[i];
    }

    arcs[middle_arc->sister].r_cap += bottle_neck;
    middle_arc->r_cap -= bottle_neck;

    // the source tree
    i = arcs[middle_arc->sister].head;
    while (true) {
        aid = parent[i];
        if (aid == TERMINAL) {
            break;
        }
        Arc& arc = arcs[aid];
        arc.r_cap += bottle_neck;
        Arc& sister = arcs[arc.sister];
        sister.r_cap -= bottle_neck;
        if (sister.r_cap == 0.0) {
            // set_orphan_front(i);
        }
    }
    tr_cap[i] -= bottle_neck;
    if (tr_cap[i] == 0.0) {
        // set_orphan_front(i);
    }

    // the sink tree
    i = middle_arc->head;
    while (true) {
        aid = parent[i];
        if (aid == TERMINAL) {
            break;
        }
        Arc& arc = arcs[aid];
        arc.r_cap -= bottle_neck;
        Arc& sister = arcs[arc.sister];
        sister.r_cap += bottle_neck;
        if (arc.r_cap == 0.0) {
            // set_orphan_front(i);
        }
    }
    tr_cap[i] += bottle_neck;
    if (tr_cap[i] == 0.0) {
        // set_orphan_front(i);
    }
    flow += bottle_neck;
}

double GraphCut::max_flow() {
    uint32_t current_node = INVALID;
    while (true) {
        uint32_t i = current_node;
        if (i != INVALID) {
            next[i] = INVALID;
            if (parent[i] == INVALID) {
                i = INVALID;
            }
        }
        if (i == INVALID) {
            if ((i = next_active()) == INVALID) {
                break;
            }
        }
        uint32_t aid = INVALID;
        if (!is_sink[i]) {
            aid = first_arc[i];
            while (aid != INVALID) {
                const Arc& a = arcs[aid];
                if (a.r_cap != 0.0) {
                    const uint32_t j = a.head;
                    if (parent[j] == INVALID) {
                        is_sink[j] = false;
                        parent[j] = a.sister;
                        ts[j] = ts[i];
                        set_active(j);
                    } else if (is_sink[j]) {
                        break;
                    } else if (ts[j] <= ts[i] && dist[j] > dist[i]) {
                        parent[j] = a.sister;
                        ts[j] = ts[i];
                        dist[j] = dist[i] + 1;
                    }
                }
                aid = a.next;
            }
        } else {
            aid = first_arc[i];
            while (aid != INVALID) {
                const Arc& a = arcs[aid];
                if (a.r_cap != 0.0) {
                    const uint32_t j = a.head;
                    if (parent[j] == INVALID) {
                        is_sink[j] = true;
                        parent[j] = a.sister;
                        ts[j] = ts[i];
                        set_active(j);
                    } else if (!is_sink[j]) {
                        break;
                    } else if (ts[j] <= ts[i] && dist[j] > dist[i]) {
                        parent[j] = a.sister;
                        ts[j] = ts[i];
                        dist[j] = dist[i] + 1;
                    }
                }
                aid = a.next;
            }
        }
        time += 1;
        if (aid != INVALID) {
            Arc& a = arcs[aid];
            next[i] = i;
            current_node = i;
            augment(&a);
        }
    }
    return flow;
}
