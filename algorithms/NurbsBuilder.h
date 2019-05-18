#include <utility>

//
// Created by Andrey Matveev on 2019-05-17.
//

#ifndef OPENGL_NURBSBUILDER_H
#define OPENGL_NURBSBUILDER_H

#include <vector>
#include "../math/point.h"
#include "../kdtree.h"

struct Nurbs {
    std::vector<Vertex> grid;
};

class NurbsBuilder {
    std::vector<MyPoint> myPoints; // for kdtree
    std::vector<Vertex> verts;
    std::vector<int> vertsClasterId;
    std::vector<std::vector<int>> clastersVerts;
//    std::vector<Nurbs> vNurbs;
    kdt::KDTree<MyPoint> kdTree;
    std::vector<Triangle> triangles;
    std::vector<Vertex> newVertices;
    std::vector<std::vector<float>> vNurbs;

public:

    explicit NurbsBuilder(std::vector<Vertex> verts): verts(std::move(verts)) {}

    void buildKdTreeFromVerts() {
        for (auto v : verts) {
            myPoints.emplace_back(v.x, v.y, v.z);
        }
        kdTree.build(myPoints);
    }

    void clasterizePoints(int wantNeib = 50, int restPart = 10) {
        buildKdTreeFromVerts();

        vertsClasterId.assign(verts.size(), -1);

        // TODO: find nearest not used points to current (need to find 36 this points) (maybe shuffle points at first)
        for (int i = 0; i < verts.size(); ++i) {
            if (vertsClasterId[i] != -1) continue;
            auto nbs = kdTree.knnSearch(myPoints[i], wantNeib);
            std::vector<int> curClaster;
            for (int id : nbs) {
                if (vertsClasterId[i] == -1) {
                    curClaster.push_back(id);
                }
            }
            int rest = wantNeib - (int)curClaster.size();
            if (rest > wantNeib / restPart) {
                continue;
            }
            for (int j : curClaster) {
                vertsClasterId[j] = (int)clastersVerts.size();
            }
            clastersVerts.push_back(curClaster);
        }
    }

    std::tuple<int, int, int> find3Points(const std::vector<int> &claster) {
        if (claster.size() < 4) return {-1, -1, -1}; // 4 - just from ceiling
        int i1 = claster[0], i2 = claster[1];
        float maxDist = getDist(verts[i1], verts[i2]);
        for (int i = 0; i < claster.size(); ++i) {
            for (int j = i + 1; j < claster.size(); ++j) {
                float curDist = getDist(verts[claster[i]], verts[claster[j]]);
                if (curDist > maxDist) {
                    i1 = claster[i];
                    i2 = claster[j];
                    maxDist = curDist;
                }
            }
        }

        auto vm = getMid(myPoints[i1], myPoints[i2]).toVertex();
        maxDist = 0;
        int i3 = -1;
        for (int i = 0; i < claster.size(); ++i) {
            if (claster[i] == i1 || claster[i] == i2) continue;
            float curDist = getDist(verts[claster[i]], vm);
            if (i3 == -1 || curDist > maxDist) {
                i3 = claster[i];
                maxDist = curDist;
            }
        }

        return {i1, i2, i3};
    }

    void buildNurbsFromCluster(const std::vector<int> &claster) {
        auto tr = find3Points(claster);

        Vertex a = verts[std::get<0>(tr)];
        Vertex b = verts[std::get<1>(tr)];
        Vertex c = verts[std::get<2>(tr)];
        Vertex ab = (b - a);
        float len = ab.len();
        Vertex ac = (c - a);
        Vertex norm = ab ^ ac;
        ac = norm ^ ab;
        ab.normalize();
        ac.normalize();
        norm.normalize();
        triangles.emplace_back(a, b, c);

        int num = 6;
        Vertex start = a - (ac * (len / 2));
        float step = len / (num - 1);
        std::vector<float> nurb;
        for (int i = 0; i < num; ++i) {
            for (int j = 0; j < num; ++j) {
                Vertex cur = start + (ab * step) * i + (ac * step) * j;

                int ind = kdTree.nnSearch(cur.toMyPoint());
                Vertex near = verts[ind];

                float koef = (near - cur) * norm;
                cur = cur + norm * koef;

                cur = near;
                cur.r = 0.1f;
                cur.g = 1.0f;
                cur.b = 0.1f;
                newVertices.push_back(cur);
                nurb.push_back(cur.x);
                nurb.push_back(cur.y);
                nurb.push_back(cur.z);
            }
        }
        vNurbs.push_back(nurb);
    }

    void buildNurbs() {
        clasterizePoints();

        for (auto& cl : clastersVerts) {
            buildNurbsFromCluster(cl);
        }
    }

    std::vector<Triangle> getTriangles() {
        return triangles;
    }

    std::vector<Vertex> getNewVertices() {
        return newVertices;
    }

    std::vector<std::vector<float>> getNurbs() {
        return vNurbs;
    }
};


#endif //OPENGL_NURBSBUILDER_H
