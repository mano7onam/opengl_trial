//
// Created by Andrey Matveev on 2019-04-30.
//

#ifndef OPENGL_SETPAIRS_H
#define OPENGL_SETPAIRS_H

#include <set>
#include <map>
#include <algorithm>

using std::set;
using std::map;
using std::pair;

class SetPairs {
    set<pair<int, int>> S;
    map<int, int> M;

public:

    bool havePair(pair<int, int> p) {
        return static_cast<bool>(S.count(p));
    }

    bool havePair(int a) {
        return static_cast<bool>(M.count(a));
    }

    bool isEmpty() {
        bool empt = S.empty() || static_cast<int>(S.size()) == 0;
        std::cerr << S.size() << std::endl;
        return empt;
    }

    void addPair(pair<int, int> p) {
        if (!havePair(p)) {
            S.insert(p);
            S.insert({p.second, p.first});
            M[p.first]++;
            M[p.second]++;
        }
    }

    void erasePair(pair<int, int> p) {
        if (havePair(p)) {
            S.erase(p);
            S.erase(p);
            M[p.first]--;
            M[p.second]--;
        }
    }

    pair<int, int> getFirstPair() {
        return *S.begin();
    }
};

#endif //OPENGL_SETPAIRS_H
