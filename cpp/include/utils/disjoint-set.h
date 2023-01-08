#pragma once
#include <numeric>
#include <vector>

template <typename Index>
class DisjointSet {
  private:
    std::vector<Index> parent;
    std::vector<Index> rank;

  public:
    Index n_groups;

    DisjointSet(const Index n) : n_groups{n} {
        parent.resize(n);
        std::iota(parent.begin(), parent.end(), 0);
        rank.resize(n, 0);
    }

  private:
    void link(const Index x, const Index y) {
        if (x == y) {
            return;
        }

        if (rank[x] > rank[y]) {
            parent[y] = x;
        } else {
            parent[x] = y;
            if (rank[x] == rank[y]) {
                rank[y] += 1;
            }
        }
        n_groups -= 1;
    }

  public:
    void merge(const Index x, const Index y) { link(find_set(x), find_set(y)); }

    Index find_set(const Index x) {
        if (x != parent[x]) {
            parent[x] = find_set(parent[x]);
        }
        return parent[x];
    }
};
