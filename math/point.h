//
// Created by Andrey Matveev on 2019-05-17.
//
#pragma once

#ifndef OPENGL_POINT_H
#define OPENGL_POINT_H

#include <array>
#include <cmath>

class MyPoint;

struct Vertex
{
    float x, y, z, w;
    float r, g, b, a;

    Vertex(float x = 0, float y = 0, float z = 0): x(x), y(y), z(z) {}

    Vertex operator + (Vertex v) {
        return {x + v.x, y + v.y, z + v.z};
    }

    Vertex operator - (Vertex v) {
        return {x - v.x, y - v.y, z - v.z};
    }

    Vertex operator ^ (Vertex v) {
        float nx = (y * v.z) - (z * v.y);
        float ny = -((x * v.z) - (z * v.x));
        float nz = (x * v.y) - (y * v.x);
        return {nx, ny, nz};
    }

    Vertex operator * (float f) {
        return Vertex(x * f, y * f, z * f);
    }

    float operator * (Vertex v) {
        return x * v.x + y * v.y + z * v.z;
    }

    float len() {
        return sqrt(x * x + y * y + z * z);
    }

    void normalize() {
        float l = len();
        x /= l;
        y /= l;
        z /= l;
    }

    MyPoint toMyPoint();
};

class MyPoint : public std::array<float, 3>
{
public:

    // dimension of space (or "k" of k-d tree)
    // KDTree class accesses this member
    static const int DIM = 3;

    // the constructors
    MyPoint() {}
    MyPoint(float x, float y, float z)
    {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    }

    Vertex toVertex() {
        return Vertex((*this)[0], (*this)[1], (*this)[2]);
    }
};

MyPoint Vertex::toMyPoint() {
    return MyPoint(x, y, z);
}

#define sqr(x) ((x) * (x))

float getDist(Vertex a, Vertex b) {
    return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y) + sqr(a.z - b.z));
}

//float getDist(MyPoint a, MyPoint b) {
//    return sqrt(sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]));
//}

float getS(Vertex a, Vertex b, Vertex c) {
    return abs(a.x - b.x) + abs(a.y - b.y) + abs(a.z - b.z) +
           abs(a.x - c.x) + abs(a.y - c.y) + abs(a.z - c.z) +
           abs(b.x - c.x) + abs(b.y - c.y) + abs(b.z - c.z);
}

MyPoint getMid(MyPoint a, MyPoint b) {
    return MyPoint((a[0] + b[0]) * 0.5f, (a[1] + b[1]) * 0.5f, (a[2] + b[2]) * 0.5f);
}

struct Triangle {
    Vertex a, b, c;

    Triangle(Vertex a, Vertex b, Vertex c): a(a), b(b), c(c) {}
};

#endif //OPENGL_POINT_H
