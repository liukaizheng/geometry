#include <graphcut/graphcut.h>

GraphCut::GraphCut(const uint32_t n, const double* source_cap, const double* sink_cap)
    : n_nodes(n), TERMINAL(INVALID - 1), ORPHAN(INVALID) {
    dist.resize(n, 1);
    is_sink.resize(n);
    next.resize(n, INVALID);
    next_orphan.resize(n, INVALID);
    parent.resize(n, TERMINAL);
    ts.resize(n, 0);
    tr_cap.resize(n, 0.0);
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

void GraphCut::process_source_orphan(const uint32_t i) {
    uint32_t a0 = first_arc[i];
    uint32_t a0_min = INVALID;
    uint32_t d;
    uint32_t d_min = INVALID;
    for (a0 = first_arc[i]; a0 != INVALID; a0 = arcs[a0].next) {
        if (arcs[arcs[a0].sister].r_cap == 0.0) {
            continue;
        }
        uint32_t j = arcs[a0].head;
        if (is_sink[j] || parent[j] == INVALID) {
            continue;
        }
        // checking the origin of j
        d = 0;
        while (true) {
            if (ts[j] == time) {
                d += dist[j];
                break;
            }
            const uint32_t a = parent[j];
            d++;
            if (a == TERMINAL) {
                ts[j] = time;
                dist[j] = 1;
                break;
            }
            if (a == ORPHAN) {
                d = INVALID;
                break;
            }
            j = arcs[a].head;
        }
        if (d != INVALID) { // j originates from the source - done
            if (d < d_min) {
                a0_min = a0;
                d_min = d;
            }
            // set marks along the path
            for (j = arcs[a0].head; ts[j] != time; j = arcs[parent[j]].head) {
                ts[j] = time;
                dist[j] = d--;
            }
        }
    }

    parent[i] = a0_min;
    if (parent[i] != INVALID) {
        ts[i] = time;
        dist[i] = d_min + 1;
    } else {
        // process neighbors
        for (a0 = first_arc[i]; a0 != INVALID; a0 = arcs[a0].next) {
            const uint32_t j = arcs[a0].head;
            const uint32_t a = parent[j];
            if (!is_sink[j] && a != INVALID) {
                if (arcs[arcs[a0].sister].r_cap != 0.0) {
                    set_active(j);
                }
                if (a != TERMINAL && a != ORPHAN && arcs[a].head == i) {
                    set_orphan_rear(j); // add j to the end of the adoption list
                }
            }
        }
    }
}

void GraphCut::process_sink_orphan(const uint32_t i) {
    uint32_t a0 = first_arc[i];
    uint32_t a0_min = INVALID;
    uint32_t d;
    uint32_t d_min = INVALID;
    // trying to find a new parent
    for (a0 = first_arc[i]; a0 != INVALID; a0 = arcs[a0].next) {
        if (arcs[a0].r_cap == 0.0) {
            continue;
        }
        uint32_t j = arcs[a0].head;
        if (!is_sink[j] || parent[j] == INVALID) {
            continue;
        }
        // checking the origin of j
        d = 0;
        while (true) {
            if (ts[j] == time) {
                d += dist[j];
                break;
            }
            const uint32_t a = parent[j];
            d++;
            if (a == TERMINAL) {
                ts[j] = time;
                dist[j] = 1;
                break;
            }
            if (a == ORPHAN) {
                d = INVALID;
                break;
            }
            j = arcs[a].head;
        }
        if (d != INVALID) { // j originates from the sink - done
            if (d < d_min) {
                a0_min = a0;
                d_min = d;
            }
            // set marks along the path
            for (j = arcs[a0].head; ts[j] != time; j = arcs[parent[j]].head) {
                ts[j] = time;
                dist[j] = d--;
            }
        }
    }

    parent[i] = a0_min;
    if (parent[i] != INVALID) {
        ts[i] = time;
        dist[i] = d_min + 1;
    } else {
        // no parent is found, process neighbors
        for (a0 = first_arc[i]; a0 != INVALID; a0 = arcs[a0].next) {
            const uint32_t j = arcs[a0].head;
            const uint32_t a = parent[j];
            if (is_sink[j] && a != INVALID) {
                if (arcs[a0].r_cap != 0.0) {
                    set_active(j);
                }
                if (a != TERMINAL && a != ORPHAN && arcs[a].head == i) {
                    set_orphan_rear(j); // add j to the end of the adoption list
                }
            }
        }
    }
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
        i = arcs[aid].head;
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
        i = arc.head;
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
            set_orphan_front(i);
        }
        i = arc.head;
    }
    tr_cap[i] -= bottle_neck;
    if (tr_cap[i] == 0.0) {
        set_orphan_front(i);
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
            set_orphan_front(i);
        }
        i = arc.head;
    }
    tr_cap[i] += bottle_neck;
    if (tr_cap[i] == 0.0) {
        set_orphan_front(i);
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
                        dist[j] = dist[i] + 1;
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
                if (arcs[a.sister].r_cap != 0.0) {
                    const uint32_t j = a.head;
                    if (parent[j] == INVALID) {
                        is_sink[j] = true;
                        parent[j] = a.sister;
                        ts[j] = ts[i];
                        dist[j] = dist[i] + 1;
                        set_active(j);
                    } else if (!is_sink[j]) {
                        aid = a.sister;
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
            while ((i = orphan_first) != INVALID) {
                const uint32_t next_node = next_orphan[i];
                next_orphan[i] = INVALID;
                while ((i = orphan_first) != INVALID) {
                    orphan_first = next_orphan[i];
                    if (orphan_first == INVALID) {
                        orphan_last = INVALID;
                    }
                    if (is_sink[i]) {
                        process_sink_orphan(i);
                    } else {
                        process_source_orphan(i);
                    }
                }
                orphan_first = next_node;
            }
        } else {
            current_node = INVALID;
        }
    }
    return flow;
}
