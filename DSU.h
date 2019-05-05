//
// Created by Andrey Matveev on 2019-05-05.
//

#ifndef OPENGL_DSU_H
#define OPENGL_DSU_H

#include <vector>

using std::vector;

class DSU {
    vector<int> parent;
    vector<int> rank;

public:

    DSU(const vector<int> &parent) : parent(parent) {
        rank.assign(parent.size(), 0);
    }

    DSU(size_t sz) {
        rank.assign(sz, 0);
        parent.assign(sz, 0);
        for (int i = 0; i < parent.size(); ++i) {
            parent[i] = i;
        }
    }

    DSU() {}

    void init(size_t sz) {
        rank.assign(sz, 0);
        parent.assign(sz, 0);
        for (int i = 0; i < parent.size(); ++i) {
            parent[i] = i;
        }
    }

    int findSet (int v) {
        if (v == parent[v]) {
            return v;
        }
        return parent[v] = findSet (parent[v]);
    }

    void unionSets (int a, int b) {
        a = findSet (a);
        b = findSet (b);
        if (a != b) {
            if (rank[a] < rank[b]) {
                std::swap(a, b);
            }
            parent[b] = a;
            if (rank[a] == rank[b]) {
                ++rank[a];
            }
        }
    }
};

#endif //OPENGL_DSU_H
